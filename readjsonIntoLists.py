import json

def returnLists(filename):
    """
    opens file and reads json string
    """

    indices = []
    energies = []
   
    f = open(filename)
    dic = json.load(f)
    f.close()     

    for k in dic.keys():
        if k.startswith('Surface'):
            indexstr = k.strip('Surface Energy[').strip(']').split(',')
            ind = [int(x) for x in indexstr]
            indices.append(ind)        
            energies.append(dic[k])

    return indices, energies
