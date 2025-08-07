
# nf-core/sarek pipeline public datasets run

Same setup as for `nfcore_sarek_rerun` by without DeepVarinat calling. 

## Run commands
Run in screen

Setup env
```
source /vulpes/ngi/production/latest/conf/sourceme_sthlm.sh
source activate NGI
```

Run nextflow

```
nextflow run /vulpes/ngi/production/v24.07/sw/sarek/3_4_2/ -profile uppmax --project ngi2016004 -c /vulpes/ngi/production/v24.07/conf/sarek_sthlm.config -c nextflow.config -params-file params.yaml
```



