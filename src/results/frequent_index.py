import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pathlib import Path

# Official top-level package bindings verified by the notebook
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from tfitpy.datasets.go import load_go

#ppi 
from tfitpy.datasets.ppi import PPI_DATASETS
# from tfitpy.indices.ppi import shared_partners


def _process_go_chunk(
    df_chunk,
    background_gene_list,
    data_path,
    threshold):
    """
    Worker helper executed inside isolated process blocks.
    Loads GO data arrays once from disk, then loops through chunk items.
    """
    if df_chunk.empty:
        return df_chunk

    # 1. Local isolated memory load
    go_data = load_go(data_path)
    godag = go_data["godag"]
    gene2go = go_data["gene2go"]

    bg_set = set(str(g) for g in background_gene_list)

    # 2. Reconstruct flat assoc map into the exact namespace schema expected by GOATOOLS:
    # dict of: { 'BP': { gene: {go_ids} }, 'MF': { gene: {go_ids} }, 'CC': { gene: {go_ids} } }
    ns2assoc = {'BP': {}, 'MF': {}, 'CC': {}}
    
    for gene, go_ids in gene2go.items():
        for go_id in go_ids:
            if go_id in godag:
                ns = godag[go_id].namespace
                
                # Normalize namespace string names to the shorthand keys expected
                if ns == 'biological_process':
                    ns_key = 'BP'
                elif ns == 'molecular_function':
                    ns_key = 'MF'
                elif ns == 'cellular_component':
                    ns_key = 'CC'
                else:
                    continue  # Skip unmapped root tags
                
                # Build the inner dict and set association safely
                if gene not in ns2assoc[ns_key]:
                    ns2assoc[ns_key][gene] = set()
                ns2assoc[ns_key][gene].add(go_id)

    # 3. Instantiate the Official Namespace Study Handler matching notebook step 3
    goea = GOEnrichmentStudyNS(
        pop=bg_set,
        ns2assoc=ns2assoc,  # Properly structured namespace mapping
        godag=godag,
        propagate_counts=False,
        alpha=threshold,
        methods=['fdr_bh']
    )

    go_scores = []
    go_pvalues = []
    go_classes = []

    # 4. Step through chunk rows sequentially
    for item in df_chunk['common_targets']:
        if isinstance(item, str):
            targets = [g.strip() for g in item.replace(';', ',').split(',')]
        else:
            targets = [str(g) for g in item]

        valid_targets = set(g for g in targets if g in bg_set)

        # Statistical baseline safety gate
        if len(valid_targets) < 2:
            go_scores.append(0.0)
            go_pvalues.append(1.0)
            go_classes.append("")
            continue

        # 5. Run ORA quietly via prt=None as verified in step 5a of the notebook
        results_all = goea.run_study(valid_targets, prt=None)

        # Filter out non-significant rows using the p_fdr_bh attribute lookup
        significant_results = [
            res for res in results_all if getattr(res, 'p_fdr_bh', 1.0) <= threshold
        ]

        if not significant_results:
            go_scores.append(0.0)
            go_pvalues.append(1.0)
            go_classes.append("")
        else:
            # Sort records in ascending significance order (smallest p-value first)
            significant_results.sort(key=lambda x: getattr(x, 'p_fdr_bh', 1.0))

            # Extract the minimum adjusted value
            min_p = float(getattr(significant_results[0], 'p_fdr_bh'))
            
            # Guard against log10(0) bounds
            score = -np.log10(min_p) if min_p > 0 else -np.log10(1e-300)

            # Map the functional label strings from the result record properties
            class_strings = [
              f"{res.name} ({res.GO})" for res in significant_results
            ]
            classes_str = "; ".join(class_strings)

            go_scores.append(round(score, 5))
            go_pvalues.append(round(min_p, 5) if min_p > 1e-5 else min_p)
            go_classes.append(classes_str)

    # 6. Map vector results back to a structured copy
    df_out = df_chunk.copy()
    df_out['go_score'] = go_scores
    df_out['go_pvalue'] = go_pvalues
    df_out['go_classes'] = go_classes

    return df_out


def compute_GO_ORA(
    data,
    background_gene_list,
    data_path,
    batch = 4,
    threshold = 0.05):
    """
    Partitions the data matrix into balanced sub-arrays, processes them
    concurrently across process spaces to maximize core utilization, and reassembles 
    the elements to preserve matrix indexing.
    """
    path_obj = Path(data_path)
    
    # Resolve joblib core configuration values (handling negative indices like -1)
    import os
    if batch <= 0:
        n_workers = os.cpu_count() or 1
    else:
        n_workers = int(batch)
        
    if n_workers > len(data):
        n_workers = len(data)
        
    if data.empty or n_workers <= 0:
        return data

    # Partition DataFrame arrays symmetrically based on core allocation limits
    chunks = np.array_split(data, n_workers)
    print(f"Dividing pipeline work array into {len(chunks)} blocks across {n_workers} core(s).")

    processed_chunks = Parallel(n_jobs=n_workers, backend="loky")(
        delayed(_process_go_chunk)(
            df_chunk=chunk,
            background_gene_list=background_gene_list,
            data_path=path_obj,
            threshold=threshold
        ) for chunk in chunks
    )

    # Concat contiguous chunks sequentially to match structural input row alignments
    final_df = pd.concat(processed_chunks, axis=0)
    return final_df


def _process_ppi(data, ppi_name, data_path):
    """
    Processes a specific PPI source and returns JUST the new score column
    to save memory and keep structural alignment clean.
    """
    # Load the specific network for this worker
    ppi_network = PPI_DATASETS[ppi_name]["load"](data_path)

    def calculate_row_partners(row):
        item = row["common_targets"]
        # Safe handling for strings vs pre-parsed lists
        if isinstance(item, str):
            targets = [g.strip() for g in item.split(';')]
        elif isinstance(item, (list, set)):
            targets = list(item)
        else:
            return 0  # Fallback score if empty/null

        scr, all_scr = shared_partners(targets, ppi_network)
        return scr

    # Return a single Series with the original index to ensure perfect alignment
    return pd.Series(
        data.apply(calculate_row_partners, axis=1),
        name=f'ppi_shared_partners_{ppi_name}',
        index=data.index
    )


def compute_PPI_shared_partners(data, data_path):
    """
    Runs PPI processing in parallel across different PPI sources and 
    merges columns side-by-side while preserving the original DataFrame shape.
    """
    path_obj = Path(data_path)
    ppi_sources = ["hippie", "stringdb", "biogrid"]

    # Run parallel workers to get the new columns
    new_columns = Parallel(n_jobs=len(ppi_sources))(
        delayed(_process_ppi)(
            data=data,
            ppi_name=ppi,
            data_path=path_obj,
        ) for ppi in ppi_sources
    )
    
    # Horizontally glue the new columns onto the original data
    final_df = pd.concat([data] + new_columns, axis=1)
    
    return final_df