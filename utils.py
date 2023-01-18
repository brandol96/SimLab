def write_csv(out_path, file_name, **kwargs):
    import csv
    import numpy as np
    with open(f'{out_path}{file_name}.csv', 'w') as out_file:
        file_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        cols_list = []
        for key, value in kwargs.items():
            # first I have to join key and value in a single list
            # after a little check, so I can work with arrays or list freely
            ar = np.array([])
            if type(value) == list:
                cols_list.append([key]+value)
            if type(value) == type(ar):
                value = value.tolist()
                cols_list.append([key]+value)
        rows = zip(*cols_list)
        for row in rows:
            file_writer.writerow(row)
