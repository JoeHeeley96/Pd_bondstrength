def opt_check(file, string_to_search):
    with open(file, 'r') as f:
        for line in f:
            if string_to_search in line:
                return True
    return False