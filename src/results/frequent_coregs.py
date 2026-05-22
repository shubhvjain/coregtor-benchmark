import  src.results.util as ut 

from mlxtend.preprocessing import TransactionEncoder
from mlxtend.frequent_patterns import fpgrowth


def gene_frequency_patterns(results,min_support=0.01, result_folder=None,d=None):
    """"""

    transactions = results["sources"].str.split(";").tolist()
    te = TransactionEncoder()
    te_array = te.fit(transactions).transform(transactions)
    df_encoded = pd.DataFrame(te_array, columns=te.columns_)
    frequent_itemsets = fpgrowth(df_encoded, min_support=min_support, use_colnames=True)
    frequent_itemsets["length"] = frequent_itemsets["itemsets"].apply(len)
    frequent_itemsets = frequent_itemsets[frequent_itemsets['length']>1]
    frequent_itemsets["itemsets"] = frequent_itemsets["itemsets"].apply(lambda x : ";".join(x)) 
    frequent_itemsets = frequent_itemsets.rename(columns={"itemsets":"sources"})
    frequent_itemsets["dataset"] = d
    return frequent_itemsets


def coreg_frequency(exp_name, exp_data, result_data, rerun, n_jobs=1):
    """"""   
    result_folder = get_output_path() / exp_name/"_results"/ "coregulators_frequency"
    result_folder.mkdir(parents=True, exist_ok=True)

    result_file = result_folder / f"{result_data['name']}.csv"
    if result_file.exists() and not rerun:
        raise FileExistsError(result_file)
    
    ipath = get_output_path()/exp_name/"_results"/"filter_target_genes"/f"{result_data['target_gene_list']}.json"
    with open(ipath,"r")as f:
        target_input =  json.load(f)
        # print(target_input)
    frequent_genes = []

    for d in exp_data["datasets"]:
        print(d)
        fpath = get_output_path()/exp_name/d/"coregtor"/f"result_{result_data['result_name']}.csv"
        results = pd.read_csv(fpath)
        targets = target_input["list"][d]
        results = results[results['target'].isin(targets)]
        min_support = result_data.get("min_support",0.01)
        outputs = gene_frequency_patterns(results,min_support,result_folder,d)
        frequent_genes.append(outputs)

    df = pd.concat(frequent_genes)
    # print(df)
    df["support_per"] = df["support"]*100
    df.to_csv(result_file,index=False)


def find_frequent_coregs(input, env, options, args):
    """
    """
    out_path, temp_path = ut.get_exp_path(input, env)

    n_jobs = args.njobs if args.njobs is not None else -1
    rerun = args.rerun

    cluster_folder = temp_path/"clusters"
    cluster_folder.mkdir(exist_ok=True, parents=True)

    cluster_file_list = ut.get_cluster_list_result(input,options.get("result_name","NotFound"))
    all_results = [cl["id"] for cl in input["clustering"]]
    
    result_list_input = all_results
    if cluster_file_list is not None :
        result_list_input = cluster_file_list

    result_list = []  # only that are not already done
    for r in result_list_input:
        result_file_path = cluster_folder / f"freq_coreg_{r}.csv"
        if result_file_path.exists() and not rerun:
            continue
        result_list.append(r)   
    
    for r in result_list:
        print(r)
        cluster_file = cluster_folder / f"{r}.csv"
        index_file = cluster_folder / f"freq_coreg_{r}.csv"
        data = pd.read_csv(cluster_file)
        data.rename(columns={"uid":"cluster_uid"}, inplace=True)

         


