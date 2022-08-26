
def check_results(res):
    if len(res) == 0:
        raise ValueError("Did not find any features with those names, "
                         "plasae check the names you provided")
    else:
        return res

