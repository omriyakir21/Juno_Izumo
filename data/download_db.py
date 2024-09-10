import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import requests
import paths

def download_uniparc_data():
    """
    Download the UniParc data for the human taxonomy ID (9606) in FASTA format.
    
    Parameters
    ----------
    data_path : str
        Path to the data directory.
    """
    # API URL and parameters
    url = "https://rest.uniprot.org/uniparc/stream"
    params = {
        "compressed": "true",
        "format": "fasta",
        "query": "((taxonomy_id:9606))"
    }

    # Define the download path
    download_path = os.path.join(paths.data_path, 'uniparc_Humo_Sapiens_data.gzip')

    # Make the GET request
    response = requests.get(url, params=params)
    response.raise_for_status()  # Raises an HTTPError for bad responses (4xx and 5xx)

    # Save the response content to a file
    with open(download_path, 'wb') as file:
        file.write(response.content)

    print(f"Data downloaded to {download_path}")

if __name__ == "__main__":
    download_uniparc_data()