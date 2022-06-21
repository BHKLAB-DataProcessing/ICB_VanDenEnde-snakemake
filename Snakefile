from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]
data_source  = "https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_VanDenEnde-data/main/"

rule get_MultiAssayExp:
    output:
        S3.remote(prefix + filename)
    input:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript -e \
        '
        source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/get_MultiAssayExp.R");
        saveRDS(
            get_MultiAssayExp(study = "VanDenEnde", input_dir = paste0("{prefix}", "processed")), 
            "{prefix}{filename}"
        );
        '
        """

rule format_data:
    output:
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR.csv"),
        S3.remote(prefix + "processed/cased_sequenced.csv")
    input:
        S3.remote(prefix + "download/Vandenende.RData"),
        S3.remote(prefix + "download/GSE165252_series_matrix.txt")
    resources:
        mem_mb=2000
    shell:
        """
        Rscript scripts/Format_Data.R \
        {prefix}download \
        {prefix}processed \
        """

rule download_data:
    output:
        S3.remote(prefix + "download/Vandenende.RData"),
        S3.remote(prefix + "download/GSE165252_series_matrix.txt")
    resources:
        mem_mb=2000
    shell:
        """
        wget {data_source}Vandenende.RData -O {prefix}download/Vandenende.RData
        wget {data_source}GSE165252_series_matrix.txt -O {prefix}download/GSE165252_series_matrix.txt
        """ 