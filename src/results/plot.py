"""
Plot coregulatory modules 
"""
import src.results.util as ut
from joblib import load, Parallel, delayed
import pandas as pd
from pathlib import Path
import os

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def create_network_diagram(targets_list, details_df, evidence=None, settings=None,figure_path=None):
    """
    Generates a network diagram where primary target-source links are straight lines,
    and multiple PPI + TFBS evidence sources are cleanly stacked as curved arcs next to them.
    Evidence layers can be toggled on/off individually.
    """

    default_settings = {
        'target_color': '#FD9E02',
        'source_color': '#EBEBEB',
        'target_source_edge_color': '#A9A9A9',
        'ppi_hippie_color': '#E69F00',
        'ppi_stringdb_color': '#009E73',
        'ppi_biogrid_color': '#56B4E9',
        'tfbs_edge_color': '#CC79A7',
        'node_size': 600,
        'node_shape': 'o',
        'font_size': 9,
        'figsize': (20,12),
        'layout': 'kamada_kawai',
        'layout_scale': 20,
        'edge_width':1,
        'show_ppi_hippie': True,
        'show_ppi_stringdb': True,
        'show_ppi_biogrid': True,
        'show_tfbs': True
    }
    current_settings = {**default_settings, **(settings if settings is not None else {})}

    G = nx.Graph()

    filtered_df = details_df[details_df['target'].isin(targets_list)].copy()
    if filtered_df.empty:
        print(f"No data found for the given targets: {targets_list}")
        return

    primary_edges = set()
    evidence_edges = []

    for _, row in filtered_df.iterrows():
        target_node = row['target']
        cluster_uid = row['cluster_uid']
        sources = str(row['sources']).split(';') if pd.notna(row['sources']) else []

        G.add_node(target_node, type='target')

        for source_node in sources:
            source_node = source_node.strip()
            if source_node:
                G.add_node(source_node, type='source')
                G.add_edge(source_node, target_node)
                primary_edges.add(tuple(sorted((source_node, target_node))))

        if evidence and cluster_uid in evidence and evidence[cluster_uid] is not None:
            ppi_df = evidence[cluster_uid].get('ppi')
            if isinstance(ppi_df, pd.DataFrame) and not ppi_df.empty:
                for _, ppi_row in ppi_df.iterrows():
                    n1 = ppi_row['node1']
                    n2 = ppi_row['node2']
                    ev_type = str(ppi_row['evidence_type']).strip().lower()

                    if G.has_node(n1) and G.has_node(n2):
                        edge_pair = tuple(sorted((n1, n2)))
                        if 'hippie' in ev_type:
                            evidence_edges.append({'edge': edge_pair, 'type': 'hippie'})
                        elif 'stringdb' in ev_type:
                            evidence_edges.append({'edge': edge_pair, 'type': 'stringdb'})
                        elif 'biogrid' in ev_type:
                            evidence_edges.append({'edge': edge_pair, 'type': 'biogrid'})

            tfbs_df = evidence[cluster_uid].get('tfbs')
            if isinstance(tfbs_df, pd.DataFrame) and not tfbs_df.empty:
                for _, tfbs_row in tfbs_df.iterrows():
                    source_node = tfbs_row['source']
                    if G.has_node(source_node) and G.has_node(target_node):
                        edge_pair = tuple(sorted((source_node, target_node)))
                        evidence_edges.append({'edge': edge_pair, 'type': 'tfbs'})

    node_colors = [
        current_settings['target_color'] if G.nodes[node]['type'] == 'target'
        else current_settings['source_color']
        for node in G.nodes()
    ]
    node_labels = {node: node for node in G.nodes()}

    if current_settings['layout'] == 'spring':
        pos = nx.spring_layout(G, k=1.2, iterations=100, scale=current_settings['layout_scale'], seed=42)
    elif current_settings['layout'] == 'circular':
        pos = nx.circular_layout(G)
    elif current_settings['layout'] == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G, scale=current_settings['layout_scale'])
    else:
        pos = nx.fruchterman_reingold_layout(G, scale=current_settings['layout_scale'])

    plt.figure(figsize=current_settings['figsize'])

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        node_size=current_settings['node_size'],
        node_shape=current_settings['node_shape']
    )

    nx.draw_networkx_edges(
        G,
        pos,
        edgelist=list(primary_edges),
        edge_color=current_settings['target_source_edge_color'],
        width=current_settings['edge_width']
    )

    type_styles = {
        'hippie': {'color': current_settings['ppi_hippie_color'], 'style': 'dashed'},
        'stringdb': {'color': current_settings['ppi_stringdb_color'], 'style': 'dashed'},
        'biogrid': {'color': current_settings['ppi_biogrid_color'], 'style': 'dashed'},
        'tfbs': {'color': current_settings['tfbs_edge_color'], 'style': 'dotted'}
    }

    evidence_visibility = {
        'hippie': current_settings['show_ppi_hippie'],
        'stringdb': current_settings['show_ppi_stringdb'],
        'biogrid': current_settings['show_ppi_biogrid'],
        'tfbs': current_settings['show_tfbs']
    }

    seen_evidence = set()
    pair_counters = {}
    active_types = set()

    for ev in evidence_edges:
        pair = ev['edge']
        t = ev['type']

        if not evidence_visibility.get(t, True):
            continue

        edge_key = (pair, t)
        if edge_key in seen_evidence:
            continue

        seen_evidence.add(edge_key)
        active_types.add(t)

        if pair not in pair_counters:
            pair_counters[pair] = 0
        else:
            pair_counters[pair] += 1

        current_track_idx = pair_counters[pair]
        cfg = type_styles[t]
        calculated_rad = 0.10 + (current_track_idx * 0.10)

        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=[pair],
            edge_color=cfg['color'],
            width=current_settings['edge_width'] * 1.3,
            style=cfg['style'],
            connectionstyle=f"arc3,rad={calculated_rad}",
            arrows=True,
            arrowstyle='-'
        )

    nx.draw_networkx_labels(
        G,
        pos,
        labels=node_labels,
        font_size=current_settings['font_size']
    )

    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=current_settings['source_color'],
                   markersize=10, label='Source'),
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=current_settings['target_color'],
                   markersize=10, label='Target'),
        plt.Line2D([0], [0], color=current_settings['target_source_edge_color'],
                   lw=2, label='Target-Source Link')
    ]

    if 'hippie' in active_types:
        legend_elements.append(
            plt.Line2D([0], [0], color=current_settings['ppi_hippie_color'],
                       lw=2, linestyle='--', label='PPI: HIPPIE')
        )
    if 'stringdb' in active_types:
        legend_elements.append(
            plt.Line2D([0], [0], color=current_settings['ppi_stringdb_color'],
                       lw=2, linestyle='--', label='PPI: STRINGdb')
        )
    if 'biogrid' in active_types:
        legend_elements.append(
            plt.Line2D([0], [0], color=current_settings['ppi_biogrid_color'],
                       lw=2, linestyle='--', label='PPI: BioGRID')
        )
    if 'tfbs' in active_types:
        legend_elements.append(
            plt.Line2D([0], [0], color=current_settings['tfbs_edge_color'],
                       lw=2, linestyle=':', label='TFBS Binding')
        )

    plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.0, 1.0))
    plt.axis('off')

    if figure_path:
        plt.savefig(figure_path, format='pdf', bbox_inches='tight')
    # plt.show()


def generate_network_plot(input, env, options, args):
    """
    """
    out_path, temp_path = ut.get_exp_path(input, env)
    rerun = args.rerun

    # read the info
    indices_path = out_path / f"{options['cluster_name']}_indices.csv"
    evidence_path = out_path / f"{options['cluster_name']}_evidence.pkl"
    indices = pd.read_csv(indices_path)
    evidence = load(evidence_path)

    fig_path = out_path / f"{options['result_name']}.pdf"
    
    if fig_path.exists() and not rerun:
        print("file already exists")
        return
    
    targets = options["targets"]
    settings = options["settings"]
    create_network_diagram(targets, indices,evidence,settings,fig_path)

    



