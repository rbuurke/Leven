import pandas as pd


def expand(df: pd.DataFrame, output_filename: str):
    split_df_loc = df.iloc[:, 0]
    split_df_words = df.iloc[:, 1:]

    long_list = []

    for loc in enumerate(split_df_loc):
        for word_col in split_df_words:
            long_list.append([loc[1],
                              word_col,
                              split_df_words[word_col][loc[0]]])

    long_df = pd.DataFrame(long_list)
    long_df.columns = ['location', 'word', 'transcription']

    long_df.to_csv(output_filename, sep='\t', index=False)

    return 0
