import time
import pandas as pd
import requests as api
from Bio import Entrez, SeqIO
import logging
import ssl

ncbi_api_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
Entrez.email = "derektshiabuila@gmail.com"
Entrez.api_key = "5015c5c694ad358d33e88cd866358909f508"

# NCBI RESPONSE DATA SCHEMES
from typing import List, Optional
from pydantic import BaseModel

#########################
# DTO SUBMODELS
#########################


class VirusIsolate(BaseModel):
    name: Optional[str] = None
    source: Optional[str] = None
    collection_date: Optional[str] = None


class Lineage(BaseModel):
    tax_id: Optional[int] = None
    name: Optional[str] = None


class VirusHost(BaseModel):
    tax_id: Optional[int] = None
    sci_name: Optional[str] = None
    organism_name: Optional[str] = None
    common_name: Optional[str] = None
    lineage: List[Lineage] = None
    strain: Optional[str] = None
    pangolin_classification: Optional[str] = None



class Virus(BaseModel):
    tax_id: Optional[int] = None
    sci_name: Optional[str] = None
    organism_name: Optional[str] = None
    common_name: Optional[str] = None
    lineage: Optional[List[Lineage]] = None
    strain: Optional[str] = None
    pangolin_classification: Optional[str] = None


class ReportLocation(BaseModel):
    geographic_location: Optional[str] = None
    geographic_region: Optional[str] = None
    usa_state: Optional[str] = None


class ReportNucleotide(BaseModel):
    seq_id: Optional[str] = None
    accession_version: Optional[str] = None
    title: Optional[str] = None
    sequence_hash: Optional[str] = None


class ReportSubmitter(BaseModel):
    names: Optional[List[str]] = None
    affiliation: Optional[str] = None
    country: Optional[str] = None


#########################
# DTO MODELS
#########################

class SequenceLengthRange(BaseModel):
    min: Optional[int] = None
    max: Optional[int] = None

class VirusFilter(BaseModel):
    accessions: Optional[List[str]] = None
    taxon: Optional[str] = None
    refseq_only: Optional[bool] = None
    annotated_only: Optional[bool] = None
    released_since: Optional[str] = None             # should be formatted as '2020-04-01T00:00:00.000Z'
    updated_since: Optional[str] = None              # should be formatted as '2020-04-01T00:00:00.000Z'
    host: Optional[str] = None
    pangolin_classification: Optional[str] = None
    geo_location: Optional[str] = None
    usa_state: Optional[str] = None
    complete_only: Optional[bool] = None
    min: Optional[int] = None
    max: Optional[int] = None



class VirusReport(BaseModel):
    accession: Optional[str] = None
    is_complete: Optional[bool] = None
    is_annotated: Optional[bool] = None
    lab_host: Optional[str] = None
    is_lab_host: Optional[bool] = None
    is_vaccine_strain: Optional[bool] = None
    source_database: Optional[str] = None
    protein_count: Optional[int] = None
    update_date: Optional[str] = None
    release_date: Optional[str] = None
    nucleotide_completeness: Optional[str] = None
    completeness: Optional[str] = None
    length: Optional[int] = None
    gene_count: Optional[int] = None
    mature_peptide_count: Optional[int] = None
    biosample: Optional[str] = None
    mol_type: Optional[str] = None
    purpose_of_sampling: Optional[str] = None
    isolate: Optional[VirusIsolate] = None
    host: Optional[VirusHost] = None
    virus: Optional[Virus] = None
    bioprojects: List[str] = None
    location: Optional[ReportLocation] = None
    nucleotide: Optional[ReportNucleotide] = None
    sra_accessions: Optional[List[str]] = None
    submitter: Optional[ReportSubmitter] = None



#########################
# HTTP RESPONSES MODELS
#########################


class VirusMetadataResponse(BaseModel):
    reports: Optional[List[VirusReport]] = None
    total_count: Optional[int] = None
    next_page_token: Optional[str] = None

class VirusRequestPayload(BaseModel):
    filter: Optional[VirusFilter] = None
    page_size: Optional[int] = None
    page_token: Optional[str] = None

# Fetches the sequence data from Entrez/Genbank given a list of accession ids

ssl._create_default_https_context = ssl._create_unverified_context

def fetch_sequences_batch(accession_list:list[str], batch_size=200, output_file="/content/", sleep_time:int = 0.1) -> str:
    with open(output_file, "w") as out_handle:
        for i in range(0, len(accession_list), batch_size):
            batch_ids = accession_list[i:i + batch_size]
            try:
                handle = Entrez.efetch(db="nucleotide", id=",".join(batch_ids), rettype="fasta", retmode="text")
                seq_records = SeqIO.parse(handle, "fasta")

                # Write batch to output file
                SeqIO.write(seq_records, out_handle, "fasta")
                handle.close()
                time.sleep(sleep_time)

            except Exception as e:
                print(f"Error fetching batch {batch_ids}: {e}")

    print(f"Sequences saved to {output_file}")
    return output_file

# Fetches metadata from ncbi datasets api
def fetch_virus_metadata(payload: VirusRequestPayload = None) -> VirusMetadataResponse:

        payload = payload.model_dump(exclude_none=True, exclude_unset=True)

        response = api.post(url=f"{ncbi_api_url}/virus", json=payload).json()

        if not response:
            return []

        return VirusMetadataResponse(**response)



# Fetch all metadata from nbci database
def fetch_all_virus_metadata(filters: VirusFilter = None, start_page: str = None, page_size:int = 1000) -> List[VirusReport]:
        reports: List[VirusReport] = []

        next_page = True
        next_page_token = start_page if start_page else None

        print("Feetching NCBI metadata")

        count = 0

        print("Fetching sequences from NCBI...")

        while next_page:
            try:
                payload = VirusRequestPayload(filter=filters, page_token=next_page_token, page_size=page_size)

                response = fetch_virus_metadata(payload)

                if not response or not response.reports:
                    next_page = False

                reports.extend(response.reports)

                if not response.next_page_token:
                    next_page = False

                next_page_token = response.next_page_token

                count = 0

                time.sleep(0.1)

            except Exception as e:
                date = time.asctime(time.localtime())
                log_message = f"Error fetching the ncbi virus metadata from page {next_page_token}\n"

                print(e)
                print(count)
                count=count+1
                if count == 10:
                    raise(e)
                time.sleep(2)
                # Don't do nothing else other than printing the error, so the request will be made again

        print("Finalised!")
        return reports

# Parsing ncbi nested metadata response to unested list of python dictionaries
def parse_response_to_dictionaries(records: List[VirusReport]) -> list[dict]:

  parsed_respnse = [
      dict(
          accession_id=record.accession,
          is_complete=record.is_complete,
          is_annotated=record.is_annotated,
          isolate_name=record.isolate.name if record.isolate else None,
          isolate_source=record.isolate.source if record.isolate else None,
          isolate_collection_date=record.isolate.collection_date if record.isolate else None,
          source_database=record.source_database,
          host_tax=record.host.tax_id if record.host else None,
          host_sci_name=record.host.sci_name if record.host else None,
          host_name=record.host.organism_name if record.host else None,
          virus_tax=record.virus.tax_id if record.virus else None,
          virus_sci_name=record.virus.sci_name if record.virus else None,
          virus_name=record.virus.organism_name if record.virus else None,
          virus_lineage=record.virus.lineage[-1].name if record.virus else None,
          location=record.location.geographic_location if record.location else None,
          region=record.location.geographic_region if record.location else None,
          length=record.length,
          gene_count=record.gene_count,
          submitter_affiliation=record.submitter.affiliation if record.submitter else None,
          submitter_country=record.submitter.country if record.submitter else None,
          update_date=record.update_date,
          release_date=record.release_date
      ) for record in records
  ]

  return parsed_respnse

#HBV Pipeline

hbv_filters = VirusFilter(
    taxon="10407"
)

hbv_resutls = fetch_all_virus_metadata(filters = hbv_filters)

hbv_df = pd.DataFrame(parse_response_to_dictionaries(hbv_resutls))

#
#
#
## CHANGE THE OUTPUT PATH
hbv_df.to_csv("/Users/derektshiabuila/Documents/PHD_2024/Hepatitis_Dash/Re_ Hepatitis Virus Dashboard_vagner/Hepatitis_dash/data_auto_get/hbv_metadata.csv")


filtered_hbv_metadata = hbv_df.dropna(subset=["location", "isolate_collection_date"], axis=0)

filtered_hbv_metadata = filtered_hbv_metadata[ filtered_hbv_metadata["length"].astype(int) >= 2900 ]

hbv_ids = list(filtered_hbv_metadata["accession_id"])

print(f"HBV number of selected sequences: {len(hbv_ids)}")

#
#
#
## CHANGE THE OUTPUT PATH
fetch_sequences_batch(accession_list=hbv_ids, batch_size=50, output_file="/Users/derektshiabuila/Documents/PHD_2024/Hepatitis_Dash/Re_ Hepatitis Virus Dashboard_vagner/Hepatitis_dash/data_auto_get/hbv_sequences.fasta")

# Defining Filters
hcv_filters = VirusFilter(
    taxon="3052230",
    sequence_length_range=SequenceLengthRange(min=8000)

)

# HCV Pipeline

hcv_resutls = fetch_all_virus_metadata(filters = hcv_filters)   # Fetch ncbi metadata

hcv_df = pd.DataFrame(parse_response_to_dictionaries(hcv_resutls))    # Creates dataframe with results

#
#
#
## CHANGE THE OUTPUT PATH
hcv_df.to_csv("/Users/derektshiabuila/Documents/PHD_2024/Hepatitis_Dash/Re_ Hepatitis Virus Dashboard_vagner/Hepatitis_dash/data_auto_get/hcv_metadata.csv")                            # Export results as csv


filtered_hcv_metadata = hcv_df.dropna(subset=["location", "isolate_collection_date"], axis=0)     # Drop sequences with None countries and collection_dates

filtered_hcv_metadata = filtered_hcv_metadata[ filtered_hcv_metadata["length"].astype(int) >= 8000 ]

hcv_ids = list(filtered_hcv_metadata["accession_id"])                                              # select the outstanding sequence ids

print(f"HCV number of selected sequences: {len(hcv_ids)}")


#
#
#
## CHANGE THE OUTPUT PATH
fetch_sequences_batch(accession_list=hcv_ids, batch_size=50, output_file="/Users/derektshiabuila/Documents/PHD_2024/Hepatitis_Dash/Re_ Hepatitis Virus Dashboard_vagner/Hepatitis_dash/data_auto_get/hcv_sequences.fasta")      # fetches sequence data for the outstanding ids