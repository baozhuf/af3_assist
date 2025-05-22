import json, os, glob, argparse
import pandas as pd
import datetime as dt

def get_args():
    parser = argparse.ArgumentParser("Pipeline to summarize AlphaFold3 results")

    parser.add_argument( "--af3_out_dir", type=str, required=True, help="is the same as the AFOUT directory from the AlphaFold3 running step")
    parser.add_argument( "--summary_path", type=str, default="AF3_results_summary.csv", help="a csv file path to save the summarized AF3 results")

    args = parser.parse_args()
    return args

def extract_res_multimer(af3_out_dir, save_csv_path):
    print(af3_out_dir)
    js_path = os.path.join(af3_out_dir, '*/*/*summary_confidences.json')
    files = glob.glob(js_path, recursive=True)

    # an empty list to store all lists of iptm and ptm scores
    results = []

    for file in files:
        print(file)

        # an empty list to store iptm and ptm scores for each file
        result =[]

        # get the protein-protein pair name
        kmer = file.split("/")[-3].split("_")[-1].split("-")[0]
        
        sum_conf_js = file.split("/")[-1].split("_summ")[0].split("--VS--") # split example ".../orange1.1g012861m--VS--wp_012778596.1_summary_confidences.json"
        pathogen_nm = sum_conf_js[0].upper()
        host_nm = sum_conf_js[1] # lib seq names example egx35_rs04050

        # add the protein-protein pair name to the result list
        result.append(kmer)
        result.append(pathogen_nm)
        result.append(host_nm)

        # open file and load as json ojbect
        with open(file) as f:
            data = json.load(f)

            # add iptm score
            result.append(data['iptm'])

            # add ptm score
            result.append(data['ptm'])

            # add the ranking score (0.8iptm+0.2ptm)
            result.append(round(float(data['iptm'])*0.8 + float(data['ptm'])*0.2, 3))

        # add a result(ppnm, iptm, ptm, 0.8iptm+0.2ptm) list to the results list
        results.append(result)

    # convert results to pandas DataFrame
    df = pd.DataFrame(results, columns=["Multimer", "Pathogen_protein", "Host_protein", "ipTM", "pTM", "0.8ipTM+0.2pTM"])

    # sort values by the ranking score
    df = df.sort_values(by=["ipTM", "pTM"], ascending=False)

    # save results as csv
    df.to_csv(save_csv_path, index=False)

    return df

opts = get_args()
af3_out_dir = opts.af3_out_dir
save_csv_path = opts.summary_path

extract_res_multimer(af3_out_dir, save_csv_path)
