import pandas as pd
import numpy as np
from datetime import datetime

def marker_match(reference: pd.DataFrame, matching: pd.DataFrame, method: str, d_max: int, out_path: str) -> pd.DataFrame:
    """
    Matches markers between two dataframes based on specified criteria.

    Args:
        reference: Reference dataframe.
        matching: Matching dataframe.
        method: Matching method.
        d_max: Maximum distance for matching.
        out_path: Output file path.

    Returns:
        N/A, saves three files: matched dataframe as *.csv file, a list of retainable
        markers from reference array as *_RefSelect.csv file, and a list of retainable
        probes from matching array as *_MatSelect.csv file.
    """

    start = datetime.now()

    # Define errors for function
    test_colname = ["Name", "Chr", "Position", "BAF", "LRR_mean", "LRR_sd"]
    
    if not all(col in reference.columns for col in test_colname):
        raise ValueError("Reference is missing one of the required columns: Name, Chr, Position, BAF, LRR_mean, LRR_sd.")
    
    if not all(col in matching.columns for col in test_colname):
        raise ValueError("Matching is missing one of the required columns: Name, Chr, Position, BAF, LRR_mean, LRR_sd.")
    
    if method not in ["Distance", "BAF_delta", "LRR_mean_delta", "LRR_sd_delta"]:
        raise ValueError("Method must be one of the following: Distance, BAF_delta, LRR_mean_delta, or LRR_sd_delta.")
    
    if not isinstance(d_max, (int, float)):
        raise TypeError("D_MAX must be a numeric value.")
    
    if not isinstance(out_path, str):
        raise TypeError("OutPath must be a string value containing a valid path (without file extensions).")

    # Fetch the list of available chromosomes
    chromosomes = sorted(reference['Chr'].unique())

    # Create a holding dataframe
    results = []

    # Match markers
    for chr in chromosomes:
        print(f"Processing chromosome {chr}.")
        # Filter out chromosome
        ref = reference.loc[reference['Chr'].str.contains(str(chr))]
        mat = matching.loc[matching['Chr'].str.contains(str(chr))]

        # Match manifests based on identical position, annotate distances, deltas, and name identity
        perfect_match = pd.merge(ref, mat, on='Name', how='outer', suffixes=('.Ref', '.Mat'))
        perfect_match = perfect_match.dropna()
        perfect_match['Name.Ref'] = perfect_match['Name']
        perfect_match['Name.Mat'] = perfect_match['Name']
        perfect_match.drop('Name', axis = 1, inplace = True)
        perfect_match = perfect_match.assign(
                    Distance=np.abs(perfect_match['Position.Ref'] - perfect_match['Position.Mat']),
                    LRR_sd_delta=np.abs(perfect_match['LRR_sd.Ref'] - perfect_match['LRR_sd.Mat']),
                    LRR_mean_delta=np.abs(perfect_match['LRR_mean.Ref'] - perfect_match['LRR_mean.Mat']),
                    BAF_delta=np.abs(perfect_match['BAF.Ref'] - perfect_match['BAF.Mat']),
                    Position=np.where(perfect_match['Position.Ref'] == perfect_match['Position.Mat'], 'Same', 'Different'),
                    Name=np.where(perfect_match['Name.Ref'] == perfect_match['Name.Mat'], 'Same', 'Different'),
                )
        
        perfect_match = perfect_match[['Name.Ref', 'Chr.Ref', 'Position.Ref', 'BAF.Ref', 'LRR_mean.Ref', 'LRR_sd.Ref', 'Name.Mat', 'Chr.Mat', 'Position.Mat', 'BAF.Mat', 'LRR_mean.Mat', 'LRR_sd.Mat', 'Distance', 'LRR_sd_delta', 'LRR_mean_delta', 'BAF_delta', 'Position',  'Name']]

        if not perfect_match.empty:
            # Merge PerfectMatch with master DF
            results.append(perfect_match)

            # Purge matched SNPs from Ref and Mat data
            ref = ref[~ref['Name'].isin(perfect_match['Name.Ref'])]
            mat = mat[~mat['Name'].isin(perfect_match['Name.Mat'])]

        # Pre-filter `mat` to reduce inner loop iterations
        ref_min, ref_max = ref['Position'].min(), ref['Position'].max()
        mat = mat[(mat['Position'] >= ref_min - d_max) & (mat['Position'] <= ref_max + d_max)]

        for i in range(len(ref)):
            # Pull out SNP, filter matches by position distance criteria
            ref_pos = ref.iloc[i]['Position']
            match = mat.loc[(mat['Position'] >= ref_pos - d_max) & (mat['Position'] <= ref_pos + d_max)]

            # Match manifests based on position, annotate distances, deltas, and name identity
            if not match.empty:
                match = pd.merge(ref.iloc[[i]], match, how='cross', suffixes=('.Ref', '.Mat'))
                match = match.assign(
                    Distance=np.abs(match['Position.Ref'] - match['Position.Mat']),
                    LRR_sd_delta=np.abs(match['LRR_sd.Ref'] - match['LRR_sd.Mat']),
                    LRR_mean_delta=np.abs(match['LRR_mean.Ref'] - match['LRR_mean.Mat']),
                    BAF_delta=np.abs(match['BAF.Ref'] - match['BAF.Mat']),
                    Position=np.where(match['Position.Ref'] == match['Position.Mat'], 'Same', 'Different'),
                    Name=np.where(match['Name.Ref'] == match['Name.Mat'], 'Same', 'Different'),
                )

                match = match.sort_values(by=[method]).iloc[[0]]

                # Merge Match with master DF
                results.append(match)

                # Purge matched SNPs from Ref and Mat data
                mat = mat[~mat['Name'].isin(match['Name.Mat'])]

        print(f"Completed chromosome {chr}.\n")

    # Combine all results into a single DataFrame
    final_result = pd.concat(results, ignore_index=True)
    
    # Summarize results, save, and return dataframe
    print(f"Saving overall dataframe to {out_path}.csv")
    final_result.to_csv(f"{out_path}.csv", index=False)
    final_result[['Name.Ref']].drop_duplicates().to_csv(f"{out_path}_RefSelect.csv", index=False, header=False)
    final_result[['Name.Mat']].drop_duplicates().to_csv(f"{out_path}_MatSelect.csv", index=False, header=False)

    print("Process complete.\n")

    end = datetime.now()
    ttime = round((end-start).total_seconds()/60, 6)
    print("--- %s minutes ---" % ttime)
