
def processname(name, func):
    """
    Handle the Smaple name by munging the name accorign to a predefined function.

    :param name: the name of a sample
    :param func: te processing function to handle the sample name
    :return: the modifiedsaple output
    """
    if func == "underscore":
        return name.split("_")[0]
    elif func == "underscore2":
        return  ".".join(name.split("_")[:2])
    elif func == "underscore3":
        return  ".".join(name.split("_")[:3])
    elif func == "underscore4":
        return  ".".join(name.split("_")[:4])
    elif func == "underscore5":
        return  ".".join(name.split("_")[:5])
    elif func == "dot":
        return  name.split(".")[0]
    elif func == "dot2":
        return  ".".join(name.split(".")[:2])
    elif func == "dot3":
        return  ".".join(name.split(".")[:3])
    elif func == "dot4":
        return  ".".join(name.split(".")[:4])
    elif func == "dot5":
        return  ".".join(name.split(".")[:5])
    else:
        raise ValueError("names are split by underscores or dots")
