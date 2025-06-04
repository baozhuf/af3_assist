import json, os, glob, argparse
import pandas as pd
import datetime as dt
import zipfile

def get_args():
    parser = argparse.ArgumentParser("Pipeline to summarize AlphaFold3 results from the zip files that are downloaded from alphafoldserver.com")

    parser.add_argument( "--af3_zf_dir", type=str, required=True, help="The directiory contains running result zip files downloaded from alphafoldserver.com")
    parser.add_argument( "--summary_path", type=str, default="AF3Server_results_summary.csv", help="a csv file path to save the summarized AF3 results")

    args = parser.parse_args()
    return args

def zip_to_json(zip_file_path):
    with zipfile.ZipFile(zip_file_path, 'r') as zf:
        # There are many files within a zipfile, iterate them!
        for file_name in zf.namelist():
            # get the "*summary_confidences_0.json" file
            if file_name[-10:] == "ces_0.json":
                # Open this JSON file from the archive
                with zf.open(file_name) as jsf:
                    # Read as JSON data
                    js_data = json.load(jsf)
    return js_data

def extract_zf_res(af3_zf_dir, save_csv_path):
    
    # get all .zip files' paths as a list
    zip_path = os.path.join(af3_zf_dir, 'fold*.zip')
    zip_files = glob.glob(zip_path, recursive=True)

    # an empty list to store all lists of iptm and ptm scores
    results = []

    for zf in zip_files:

        # an empty list to store iptm and ptm scores for each zip_file
        result =[]

        # get the protein-protein pair name
        fileNM_Stem = zf.split("/")[-1].split(".zip")[0]

        # add the protein-protein pair name to the result list
        result.append(fileNM_Stem)

        # open zipfile and load as json ojbect
        data = zip_to_json(zf)

        # add iptm score
        result.append(data['iptm'])

        # add ptm score
        result.append(data['ptm'])

        # add the ranking score (0.8iptm+0.2ptm)
        result.append(round(float(data['iptm'])*0.8 + float(data['ptm'])*0.2, 3))
        
        # ranking score fetch
        result.append(data['ranking_score'])

        # add a result(ppnm, iptm, ptm, 0.8iptm+0.2ptm) list to the results list
        results.append(result)

    # convert results to pandas DataFrame
    df = pd.DataFrame(results, columns=["protein_pair_names", "ipTM", "pTM", "0.8ipTM+0.2pTM", "ranking_score"])

    # sort values by the ranking score
    df = df.sort_values(by=["ipTM", "pTM"], ascending=False)

    # save results as csv
    df.to_csv(save_csv_path, index=False)

    return df

opts = get_args()
af3_zf_dir = opts.af3_zf_dir
save_csv_path = opts.summary_path

extract_zf_res(af3_zf_dir, save_csv_path)
