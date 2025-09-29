import sys
from funcs import *
import squidpy as sq
import matplotlib.colors as mcolors


def main():
    # --- Argument Check ---
    # expect 3 items in total
    #if len(sys.argv) != 3:
     #   print("Error: This script requires exactly two file path arguments.")
      #  print(f"Usage: python {sys.argv[0]} <counts_file> <meta_file>")
       # sys.exit(1) 

    # --- Get Arguments ---
    #counts_file_path = sys.argv[1]
    #meta_file_path = sys.argv[2]

    # --- Print Information ---
    print("--- Python Script Starting ---")
    #print(f"Received Counts File: {counts_file_path}")
    #print(f"Received Metadata File: {meta_file_path}")
    print("------------------------------")

    # --- Dummy Analysis ---
    # This is where your actual data analysis logic would go.
    # For this example, we'll just pretend to do some work.
    try:
        ##########################################
        ##call the definitions here (main logic)##
        ##########################################
        
        
        ### ANALYSING AIH SAMPLES ###
        
        # --- Get dataset-specific paths ---
        paths = get_paths("AIH")
        counts = paths["counts"]
        meta = paths["meta"]
        outdirectory = paths["outdir"]

        # Load AIH samples
        adata_AIH_1607 = sq.read.nanostring(
            path=".",  # current directory
            counts_file=counts + "AIH_1607_expression_mat_fixed.csv", 
            meta_file=meta + "AIH_1607_spatial_metadata_fixed.csv")

        adata_AIH_1607.obs['cell_type'] = adata_AIH_1607.obs['final_cell_types'].astype('category')

        adata_AIH_6079 = sq.read.nanostring(
            path=".",  # current directory
            counts_file=counts + "AIH_6079_expression_mat_fixed.csv",  
            meta_file=meta + "AIH_6079_spatial_metadata_fixed.csv"
        )
        adata_AIH_6079.obs['cell_type'] = adata_AIH_6079.obs['final_cell_types'].astype('category')

 
        ### Run neighbourhood enrichment for each AIH sample at different radii ###

        z1607_df_30 = neighbourhood_enrichment(adata_AIH_1607, "AIH_1607", radius=30, save_plots=True, outdir=outdirectory)
        z6079_df_30 = neighbourhood_enrichment(adata_AIH_6079, "AIH_6079", radius=30, save_plots=True, outdir=outdirectory)

        z_avg_AIH_30 = average_neighbourhood_plot(
            [z6079_df_30, z1607_df_30],
            output_path_png=outdirectory + "AIH_neighborhood_enrichment_average_radius30.png",
            output_path_csv=outdirectory + "AIH_neighborhood_enrichment_average_radius30.csv",
            title="Average Neighborhood Enrichment AIH Samples - Radius 30 microns"
        )

        z1607_df_80 = neighbourhood_enrichment(adata_AIH_1607, "AIH_1607", radius=80, save_plots=True, outdir=outdirectory)
        z6079_df_80 = neighbourhood_enrichment(adata_AIH_6079, "AIH_6079", radius=80, save_plots=True, outdir=outdirectory)

        z_avg_AIH_80 = average_neighbourhood_plot(
            [z6079_df_80, z1607_df_80],
            output_path_png=outdirectory + "AIH_neighborhood_enrichment_average_radius80.png",
            output_path_csv=outdirectory + "AIH_neighborhood_enrichment_average_radius80.csv",
            title="Average Neighborhood Enrichment AIH Samples - Radius 80 microns"
        )

        # Run neighbourhood enrichment for radius 50 microns last as the current 
        # neighbourhood enrichment function modifies the adata object when run

        z1607_df_50 = neighbourhood_enrichment(adata_AIH_1607, "AIH_1607", radius=50, save_plots=True, outdir=outdirectory)
        z6079_df_50 = neighbourhood_enrichment(adata_AIH_6079, "AIH_6079", radius=50, save_plots=True, outdir=outdirectory)

        z_avg_AIH_50 = average_neighbourhood_plot(
            [z1607_df_50, z6079_df_50],
            output_path_png=outdirectory + "AIH_neighborhood_enrichment_average_radius50.png",
            output_path_csv=outdirectory + "AIH_neighborhood_enrichment_average_radius50.csv",
            title="Average Neighborhood Enrichment - AIH Samples"
        )
        
        # define a consistent color palette for cell types
        cell_types = [
            "B cell", "Epithelial", "Hepatic stellate", "Hepatocyte", 
            "Macrophage", "Monocyte", "NK cell", "Neutrophil", 
            "Plasma cell", "T cell", "Type 1 LSEC", "Type 2 LSEC"
        ]
        
        # Custom hex palette
        palette = [
            "#FFA6C9",  # B cell
            "#D5E8D4",  # Epithelial
            "#FFD1DC",  # Hepatic stellate
            "#AEC6CF",  # Hepatocyte
            "#F7B2AB",  # Macrophage
            "#FDFD96",  # Monocyte
            "#FFDAC1",  # NK cell
            "#99EE99",  # Neutrophil
            "#A9E2F3",  # Plasma cell
            "#B39EB5",  # T cell
            "#F4BBFF",  # Type 1 LSEC
            "#E6E6FA"   # Type 2 LSEC
        ]

        celltype_color_map = dict(zip(cell_types, palette))

        ### Get Ripley statistics for AIH samples ###
        
        sq.gr.ripley(adata_AIH_1607, cluster_key='cell_type', mode='G', max_dist=100)
        # Ensure consistent category order
        # Make sure cell_type is categorical
        adata_AIH_1607.obs["cell_type"] = adata_AIH_1607.obs["cell_type"].astype("category")

        categories = adata_AIH_1607.obs['cell_type'].cat.categories
        
        # Map colors in the same order
        palette_ordered = [celltype_color_map[ct] for ct in categories]

        cmap = mcolors.ListedColormap(palette_ordered)

        # Now pass it as a list
        sq.pl.ripley(
            adata_AIH_1607,
            cluster_key='cell_type',
            mode='G',
            palette=cmap, 
            plot_sims=True,
            save=outdirectory + "AIH_1607_ripley_G_sq.png"
        )

        sq.gr.ripley(adata_AIH_6079, cluster_key='cell_type', mode='G', max_dist=100)
        # Make sure cell_type is categorical
        adata_AIH_6079.obs["cell_type"] = adata_AIH_6079.obs["cell_type"].astype("category")

        categories = adata_AIH_6079.obs['cell_type'].cat.categories
        palette_ordered = [celltype_color_map[ct] for ct in categories]
        cmap = mcolors.ListedColormap(palette_ordered)

        sq.pl.ripley(
            adata_AIH_6079,
            cluster_key='cell_type',
            mode='G',
            palette=cmap,
            plot_sims=True,
            save=outdirectory + "AIH_6079_ripley_G_sq.png"
        )

        # Computing and plotting the average
        ripleyG_1607 = adata_AIH_1607.uns["cell_type_ripley_G"]
        ripleyG_6079 = adata_AIH_6079.uns["cell_type_ripley_G"]

        sims_1607G = process_sims_stat(ripleyG_1607, '1607')
        sims_6079G = process_sims_stat(ripleyG_6079, '6079')

        sims_allG = pd.concat([sims_1607G, sims_6079G], axis=0)

        plt.figure(figsize=(12, 8))

        df_1607G = ripleyG_1607['G_stat']
        df_6079G = ripleyG_6079['G_stat']

        df_1607G['source'] = '1607'
        df_6079G['source'] = '6079'

        df_allG = pd.concat([df_1607G, df_6079G], axis=0)

        df_summaryG = df_allG.groupby(['bins', 'cell_type']).agg(
            mean_stats=('stats', 'mean'),
            lower_stats=('stats', lambda x: np.percentile(x, 5)),
            upper_stats=('stats', lambda x: np.percentile(x, 95))
        ).reset_index()

        sims_summaryG = sims_allG.groupby(['bins']).agg(
            lower_sims=('stats', lambda x: np.percentile(x, 5)),
            upper_sims=('stats', lambda x: np.percentile(x, 95)),
            mean_sims=('stats', 'mean')
        ).reset_index()
        
        for cell_type in df_summaryG['cell_type'].unique():
            subset = df_summaryG[df_summaryG['cell_type'] == cell_type]
            color = celltype_color_map.get(cell_type, (0.5, 0.5, 0.5))
            
            plt.plot(subset['bins'], subset['mean_stats'], label=cell_type, color=color)
            plt.fill_between(subset['bins'], 
                     np.clip(subset['lower_stats'], 0, 1), 
                     np.clip(subset['upper_stats'], 0, 1), 
                     color=color, alpha=0.3)
        
        
        # Plot the average random expectation line (same for all cell types)
        plt.plot(sims_summaryG['bins'], sims_summaryG['mean_sims'], color='black', linestyle='--', alpha=0.5, label='Random expectation')

        plt.xlabel('Distance bins', fontsize=16, weight="bold")
        plt.ylabel('Average Ripley G statistic', fontsize=16, weight="bold")
        plt.title('Average Ripley G function across AIH slides with Random Expectation', fontsize=18, weight="bold")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "AIH_average_ripley_G_function_with_random_expectation.png")
        plt.close()

        # now for selected cell types (immune cells)
        # Define the subset of cell types you want to plot
        selected_cell_types = ['B cell', 'T cell', 'Macrophage', 'Monocyte', 'NK cell', 'Neutrophil', 'Plasma cell']
        
        plt.figure(figsize=(12, 8))
        
        for cell_type in selected_cell_types:
            subset = df_summaryG[df_summaryG['cell_type'] == cell_type]
            color = celltype_color_map.get(cell_type, (0.5, 0.5, 0.5))
            
            plt.plot(subset['bins'], subset['mean_stats'], label=cell_type, color=color)
            plt.fill_between(subset['bins'], 
                     np.clip(subset['lower_stats'], 0, 1), 
                     np.clip(subset['upper_stats'], 0, 1), 
                     color=color, alpha=0.3)
        
        # Plot the average random expectation line (same for all cell types)
        plt.plot(sims_summaryG['bins'], sims_summaryG['mean_sims'], color='black', linestyle='--', alpha=0.5, label='Random expectation')

        plt.xlabel('Distance bins', fontsize=16, weight="bold")
        plt.ylabel('Average Ripley G statistic', fontsize=16, weight="bold")
        plt.title('Average Ripley G function for selected cell types with Random Expectation', fontsize=18, weight="bold")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "AIH_average_ripley_G_immunecells.png")
        plt.close()
        
        # P vs NP G statistic
        adata_AIH_1607.obs['P_VS_NP'] = adata_AIH_1607.obs["P_VS_NP"].astype("category")
        sq.gr.ripley(adata_AIH_1607, cluster_key="P_VS_NP", mode="G", n_neigh=5, n_simulations=100, max_dist=100, n_steps=50, copy=False)
        sq.pl.ripley(adata_AIH_1607, cluster_key="P_VS_NP", mode ='G', plot_sims=True, save = outdirectory + "AIH_1607_ripley_G_sq_PvsNP.png")

        adata_AIH_6079.obs['P_VS_NP'] = adata_AIH_6079.obs["P_VS_NP"].astype("category")
        sq.gr.ripley(adata_AIH_6079, cluster_key="P_VS_NP", mode="G", n_neigh=5, n_simulations=100, max_dist=100, n_steps=50, copy=False)
        sq.pl.ripley(adata_AIH_6079, cluster_key="P_VS_NP", mode ='G', plot_sims=True, save = outdirectory + "AIH_6079_ripley_G_sq_PvsNP.png")

        cell_types_P_NP = ['P', 'NP']
        
        palette_2= [(0.9678, 0.4413, 0.5358),(0.9036, 0.5120, 0.1959)]
        
        celltype_color_map_P_NP = dict(zip(cell_types_P_NP, palette_2))

        ripleyG_1607_PNP = adata_AIH_1607.uns["P_VS_NP_ripley_G"]
        ripleyG_6079_PNP = adata_AIH_6079.uns["P_VS_NP_ripley_G"]

        plt.figure(figsize=(12, 8))

        df_1607G = ripleyG_1607_PNP['G_stat']
        df_6079G = ripleyG_6079_PNP['G_stat']

        df_1607G['source'] = '1607'
        df_6079G['source'] = '6079'

        df_allG = pd.concat([df_1607G, df_6079G], axis=0)

        df_summaryG = df_allG.groupby(['bins', 'P_VS_NP']).agg(
            mean_stats=('stats', 'mean'),
            lower_stats=('stats', lambda x: np.percentile(x, 5)),
            upper_stats=('stats', lambda x: np.percentile(x, 95))
        ).reset_index()

        sims_summaryG = sims_allG.groupby(['bins']).agg(
            lower_sims=('stats', lambda x: np.percentile(x, 5)),
            upper_sims=('stats', lambda x: np.percentile(x, 95)),
            mean_sims=('stats', 'mean')
        ).reset_index()
        
        for cell_type in df_summaryG['P_VS_NP'].unique():
            subset = df_summaryG[df_summaryG['P_VS_NP'] == cell_type]
            color = celltype_color_map_P_NP.get(cell_type, (0.5, 0.5, 0.5))
            
            plt.plot(subset['bins'], subset['mean_stats'], label=cell_type, color=color)
            plt.fill_between(subset['bins'], 
                     np.clip(subset['lower_stats'], 0, 1), 
                     np.clip(subset['upper_stats'], 0, 1), 
                     color=color, alpha=0.3)
        
        # Plot the average random expectation line (same for all cell types)
        plt.plot(sims_summaryG['bins'], sims_summaryG['mean_sims'], color='black', linestyle='--', alpha=0.5, label='Random expectation')
        plt.xlabel('Distance bins', fontsize=16, weight="bold")
        plt.ylabel('Average Ripley G statistic', fontsize=16, weight="bold")
        plt.title('Average Ripley G function across D slides with Random Expectation', fontsize=18, weight="bold")
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "AIH_average_ripley_G_function_with_random_expectation_P_VS_NP.png")
        plt.close()



        ### Now run neighbourhood enrichment per fov ###

        z_1607_per_fov = run_neighbourhood_enrichment_per_fov(adata_AIH_1607, "AIH_1607")
        z_6079_per_fov = run_neighbourhood_enrichment_per_fov(adata_AIH_6079, "AIH_6079")

        # Plot boxplots for each FOV
        z_1607_per_fov['pair'] = z_1607_per_fov.apply(canonical_pair, axis=1)

        # Compute median z-score for sorting
        pair_order = (
            z_1607_per_fov.groupby('pair')['zscore']
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        
        top20_pairs = pair_order[:20]
        # Optional: drop the original columns if not needed
        # z_all = z_all.drop(columns=['cell_type_1', 'cell_type_2'])

        plt.figure(figsize=(20, 8))

        sns.boxplot(data=z_1607_per_fov[z_1607_per_fov['pair'].isin(top20_pairs)], 
            x="pair", 
            y="zscore", 
            hue="sample", 
            palette=["#66C2A5"], 
            order=top20_pairs)

        plt.xticks(rotation=90, fontsize=18)  
        plt.yticks(fontsize=18)               
        plt.xlabel("Cell-Type Pairs", fontsize=18, weight="bold")
        plt.ylabel("Z-Score", fontsize=18, weight="bold")
        plt.title("Top 20 Neighborhood Enrichment Z-Scores per Cell-Type Pair in AIH 1607 (Ranked by Median)", 
                fontsize=20, weight="bold")
        plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=15, title_fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "top20_AIH_1607_neighborhood_enrichment_zscores_per_pair.png", dpi=300)
        print('top 20 figure saved')
        plt.close()

        z_6079_per_fov['pair'] = z_6079_per_fov.apply(canonical_pair, axis=1)

        pair_order_6079 = (
            z_6079_per_fov.groupby('pair')['zscore']
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        
        top20_pairs = pair_order_6079[:20]
        plt.figure(figsize=(20, 8))
        sns.boxplot(data=z_6079_per_fov[z_6079_per_fov['pair'].isin(top20_pairs)],
            x="pair", 
            y="zscore", 
            hue="sample", 
            palette=["#FC8D62"], 
            order=top20_pairs)

        plt.xticks(rotation=90, fontsize=18)  
        plt.yticks(fontsize=18)                
        plt.xlabel("Cell-Type Pairs", fontsize=18, weight="bold")
        plt.ylabel("Z-Score", fontsize=18, weight="bold")
        plt.title("Top 20 Neighborhood Enrichment Z-Scores per Cell-Type Pair in AIH 6079 (Ranked by Median)", 
                fontsize=20, weight="bold")
        plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=15, title_fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "top20_AIH_6079_neighborhood_enrichment_zscores_per_pair.png", dpi=300)
        print('top 20 figure saved')
        plt.close()
        
        # Combine the datasets
        combined_data = pd.concat([z_1607_per_fov, z_6079_per_fov], ignore_index=True)

        # Compute median z-score across both samples for consistent ordering
        pair_order = (
            combined_data.groupby('pair')['zscore']
            .median()
            .sort_values(ascending=False)
            .index.tolist()
        )
        top20_pairs=pair_order[:20]

        # Create the joint boxplot
        plt.figure(figsize=(24, 10))

        sns.boxplot(
            data=combined_data[combined_data['pair'].isin(top20_pairs)], 
            x="pair", 
            y="zscore", 
            hue="sample", 
            palette=["#66C2A5", "#FC8D62"], 
            order=top20_pairs
        )

        plt.xticks(rotation=90, fontsize=18)  
        plt.yticks(fontsize=18)               
        plt.xlabel("Cell-Type Pairs", fontsize=18, weight="bold")
        plt.ylabel("Z-Score", fontsize=18, weight="bold")
        plt.title("Top 20 Neighborhood Enrichment Z-Scores per Cell-Type Pair: AIH Samples (Ranked by Median)",fontsize=20, weight="bold")
        plt.xlabel("Cell-Type Pairs")
        plt.ylabel("Z-Score")
        plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=15, title_fontsize=16)
        plt.tight_layout()
        plt.savefig(outdirectory + "top20_AIH_joint_neighborhood_enrichment_zscores_per_pair.png", dpi=300, bbox_inches='tight')
        plt.close()

        # Combine the datasets
        combined_data = pd.concat([z_1607_per_fov, z_6079_per_fov], ignore_index=True)

        # Save combined z-scores as CSV
        combined_data.to_csv(outdirectory + "AIH_combined_zscore_per_fov.csv", index=False)
          
        
    except FileNotFoundError as e:
        print(f"Error: A file was not found.")
        print(e)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during analysis: {e}")
        sys.exit(1)

    print("--- Python Script Finished Successfully ---")

if __name__ == "__main__":
    main()
