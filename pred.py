import requests
import time
import json

def submit_job(file_path):
    url = "https://biosig.lab.uq.edu.au/csm_peptides/api/predict"
    with open(file_path, 'rb') as file:
        data = {'pep_list': file}
        response = requests.post(url, files=data)
        if response.status_code == 200:
            return response.json()['job_id']
        else:
            raise Exception("Error submitting job: {}".format(response.content))

def get_results(job_id):
    start_time = time.time()  
    while True:
        url = f"https://biosig.lab.uq.edu.au/csm_peptides/api/predict?job_id={job_id}"
        try:
            response = requests.get(url, timeout=60)
            if response.status_code == 200:
                data = response.json()
                if 'status' in data and data['status'] == 'running':
                    print("Job still running, waiting for results...")
                    time.sleep(30)  # Wait for 30 seconds before checking again
                else:
                    return data, time.time() - start_time
            else:
                raise Exception("Error retrieving results: {}".format(response.content))
        except requests.Timeout:
            print("Request timed out. Retrying...")
        except requests.ConnectionError:
            print("Connection error. Retrying...")

# Replace with the path to your FASTA file
job_id = submit_job("C:/Users/Rundead/Novelfold/Archive/fastas/archive900.fasta")
print("Job submitted, id:", job_id)

results, elapsed_time = get_results(job_id)
print(f"Time taken to get results: {elapsed_time} seconds")

print(json.dumps(results, indent=4))  # pretty print the results
