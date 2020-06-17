from django.db import models

# Create your models here.


class V1Model(models.Model):
    Pos_v4 = models.CharField(max_length=64)
    Pos_v1 = models.CharField(max_length=64)
    Pos_v1_Info = models.CharField(max_length=32)

    class Meta:
        db_table = 'v1'


class V2Model(models.Model):
    Pos_v4 = models.CharField(max_length=64)
    Pos_v2 = models.CharField(max_length=64)
    Pos_v2_Info = models.CharField(max_length=32)

    class Meta:
        db_table = 'v2'
