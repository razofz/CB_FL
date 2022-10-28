import os


def set_project_path():
    try: 
        path = os.environ["PROJECT_PATH"]
    except KeyError: 
        path= ''
    if len(path) == 0:
        with open("../project_path") as f:
            path = f.read()[:-1]
    # print(os.getcwd())
    os.chdir(path)
    # print("> path set to:", os.getcwd())


def create_dir(path):
    try:
        create = False
        if not os.path.exists(path):
            create = True
        os.makedirs(path, exist_ok=True)
        if create:
            print("> directory created:", path)
    except Exception as e:
        raise e


def print_files(files, in_or_out):
    """Prints input/output files neatly."""
    print(f"> {in_or_out}_files: ")
    for pair in files.items():
        if type(pair[1]) == list:
            print(f"\t{pair[0]}:")
            for item in pair[1]:
                print(f"\t\t{item}")
        elif type(pair[1]) == dict:
            print(f"\t{pair[0]}:")
            for subpair in pair[1].items():
                print(f"\t\t{subpair[0]}:\t{subpair[1]}")
        else:
            print(f"\t{pair[0]}:\t{pair[1]}")
