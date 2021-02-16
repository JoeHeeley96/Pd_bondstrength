
def print_to_csv(dataframe, filename):
    with open (filename +'.csv', 'w') as f:
        print(dataframe.to_csv(sep=','), file=f)
