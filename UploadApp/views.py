from django.shortcuts import render
from .models import V1Model, V2Model
from django.core.files.storage import FileSystemStorage
from django.conf import settings
import pandas as pd
import numpy as np
import os
import sys
import pickle


def UploadView(request):

    if request.method == 'POST' and request.FILES['vcf_file']:

        vcf_file = request.FILES['vcf_file']
        fs = FileSystemStorage()
        filename = fs.save(vcf_file.name, vcf_file)

        df_mat = pd.read_csv('%s' % os.path.join(
            settings.MEDIA_ROOT, filename), compression='gzip', sep='\n', header=None)

        def mask_comment(x):
            if x[0:2] == '##':
                return False
            else:
                return True

        m = df_mat.apply(lambda x: mask_comment(x[0]), axis=1)
        data = [x for x in df_mat[m][0][1:]]
        data = [x.split('\t') for x in data]
        c = df_mat[m].values[0][0].split('\t')

        vcf_mat = pd.DataFrame(data=data, columns=c)

        vcf_mat = vcf_mat.dropna(axis=0)

        a = np.array([len(x.split(',')[0]) for x in vcf_mat['ALT']])
        b = np.array([len(x) for x in vcf_mat['REF']])
        m_alt = (a < 2)
        m_ref = (b < 2)
        m_indel = ['INDEL' not in x.split(';') for x in vcf_mat['INFO']]

        m = m_alt & m_ref & m_indel
        vcf_mat_snp = vcf_mat[m]

        vcf_mat_snp['ix'] = vcf_mat_snp.apply(
            lambda x: x['#CHROM'] + '-' + str(x['POS']), axis=1)
        vcf_mat_snp_ix = vcf_mat_snp.set_index('ix')

        refversion = request.POST.get('Refversion')

        if refversion == 'Glycine max accession Williams genome assembly v1.0':
            with open('%s' % os.path.join(settings.MEDIA_ROOT, 'v1_list.txt'), 'rb') as f:
                snp_list = pickle.load(f)
            common_ix = list(set(snp_list) & set(list(vcf_mat_snp['ix'])))
            print(len(common_ix))
            common_values = V1Model.objects.filter(Pos_v1__in=common_ix).all()
            chrom = [x.Pos_v4.split('-')[0] for x in common_values]
            pos = [x.Pos_v4.split('-')[1] for x in common_values]
            info = [x.Pos_v1_Info for x in common_values]

        elif refversion == 'Glycine max accession Williams 82 genome assembly v2.0':
            with open('%s' % os.path.join(settings.MEDIA_ROOT, 'v2_list.txt'), 'rb') as f:
                snp_list = pickle.load(f)
            common_ix = list(set(snp_list) & set(list(vcf_mat_snp['ix'])))
            common_values = V2Model.objects.filter(Pos_v2__in=common_ix).all()
            chrom = [x.Pos_v4.split('-')[0] for x in common_values]
            pos = [x.Pos_v4.split('-')[1] for x in common_values]
            info = [x.Pos_v2_Info for x in common_values]

        vcf_mat_snp_ix_info = vcf_mat_snp_ix[vcf_mat_snp_ix.columns[0:9]]
        vcf_mat_snp_ix_samples = vcf_mat_snp_ix[vcf_mat_snp_ix.columns[9:]]

        vcf_mat_common = vcf_mat_snp_ix_info.loc[common_ix]

        vcf_mat_common['#CHROM'] = chrom
        vcf_mat_common['POS'] = pos
        vcf_mat_common['info'] = info

        def ReversePosUpdated(x):
            dic_cg = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}

            if x['info'] == False:

                FixREF.append(dic_cg[x['REF']])
                FixALT.append(','.join([dic_cg[i]
                                        for i in x['ALT'].split(',')]))

            else:
                FixREF.append(x['REF'])
                FixALT.append(x['ALT'])

        global FixREF, FixALT
        FixREF, FixALT = [], []

        progress = vcf_mat_common.apply(lambda x: ReversePosUpdated(x), axis=1)

        vcf_mat_common['REF'] = FixREF
        vcf_mat_common['ATL'] = FixALT

        df_result = pd.merge(vcf_mat_common[vcf_mat_common.columns[:-1]],
                             vcf_mat_snp_ix_samples.loc[common_ix], left_index=True, right_index=True, how='left')

        df_result['POS'] = df_result['POS'].apply(lambda x: int(x))
        df_result.sort_values(by=['#CHROM', 'POS'], inplace=True)
        m = df_result.apply(lambda x: x['#CHROM'].split('_')[
                            0] == 'scaffold', axis=1)
        df_scafold = df_result[m]
        df_scafold['int'] = df_scafold.apply(
            lambda x: int(x['#CHROM'].split('_')[1]), axis=1)
        df_scafold.sort_values(by=['int', 'POS'], inplace=True)
        df = pd.concat([df_result[~m], df_scafold[df_scafold.columns[:-1]]])

        header = '''##reference=glyma.Wm82.gnm4.4PTR.genome_main.fna\n'''
        output_VCF = '%s' % os.path.join(settings.MEDIA_ROOT, 'output.vcf')

        with open(output_VCF, 'w') as vcf:
            vcf.write(header)

        df.to_csv(output_VCF, sep='\t', mode='a', index=False)

        uploaded_file_url = fs.url('output.vcf')

        fs.delete(vcf_file.name)

        obj_output = {
            'uploaded_file_url': uploaded_file_url,
        }

        return render(request, 'upload.html', obj_output)

    return render(request, 'upload.html')
