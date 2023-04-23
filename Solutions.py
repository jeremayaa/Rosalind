
def DNA(string=None, filename=None, save=None):
    """
    Count the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in string
    Parameters:
    ----------
    filename : str, optional
        The name of the file to load data from.
    string : str, optional
        A string containing the data to be loaded.
    save_filename : str, optional
        The name of the file to save the loaded data to. If provided, the data will be saved to a file with
        this name. Default is None.

    Returns:
    -------
    str
        Four integers 
    """
    if filename is not None:
        with open(filename, 'r') as file:
            string = file.read()
    
    elif string is not None:
        pass

    else: 
        raise ValueError('Either filename or string must be provided')
    
    A = string.count('A')
    C = string.count('C')
    G = string.count('G')
    T = string.count('T')

    if save is not None:
        with open(save, 'w') as file:
            file.write(f'{A} {C} {G} {T}')

    return(A, C, G, T)


def RNA(string=None, filename=None, save=None):

    if filename is not None:
        with open(filename, 'r') as file:
            string = file.read()
    elif string is not None:
        pass
    else:
        raise ValueError('Either filename or string must be provided')
    
    string = string.replace('T', 'U')

    if save is not None:
        with open(save, 'w') as file:
            file.write(string)

    return string


def REVC(string=None, filename=None, save=None):
    """
    Returns reverse compliment of DNA string
    """

    if filename is not None:
        with open(filename, 'r') as file:
            string = file.read()

    elif string is not None:
        pass

    else:
        raise ValueError('Eihter string or filename must be provided')


    # Main algorithm
    output = ''
    for i in range(0, len(string)):
        if string[i] == 'T':
            output+='A'
        if string[i] == 'A':
             output+='T'
        if string[i] == 'G':
             output+='C'
        if string[i] == 'C':
             output+='G'

    output=output[::-1]

    if save is not None:
        with open(save, 'w') as file:
            file.write(output)

    return output


def FIB(n=None, k=None):
    """
    Fibonacci sequence. Returns number of pairs after n months, i every pair produces k-offspirng
    """
    a = 1
    b = 0
    for _ in range(0, n):
        current = a+b
        a = k*b
        b = current
        current = 0
    return b


def Fasta_read(filename):
    sequences = []
    headers = []
    with open(filename, "r") as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                headers.append(line.strip().lstrip('>'))
                sequences.append(seq)
                seq = ""
            else:
                seq += line.strip()
        sequences.append(seq)
    del sequences[0]

    return sequences, headers


def GC(filename):
    sequences, headers = Fasta_read(filename)
    current_MAX = 0
    number = 0
    for s in range(0, len(sequences)):
        current_GC = 0
        for i in sequences[s]:
            if i=='G' or i=='C':
                current_GC += 1 
        current_GC = current_GC/len(sequences[s])

        if current_GC>current_MAX:
            current_MAX=current_GC
            number = s

    return headers[number], current_MAX*100

    
def HAMM(seq1=None, seq2=None, filename=None):
    if filename is not None:
        with open(filename, 'r') as file:
            seq1, seq2 = file.read().split('\n')
    elif seq1 is not None and seq2 is not None:
        pass
    else:
        raise ValueError('provide both sequences or valid filename')

    mutations = 0
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            mutations+=1
        
    return mutations


def IPRB(k, m, n):
    T = k+m+n
    return 1-(m*n + 0.25*m*(m-1) + n*(n-1))/(T*(T-1))


import pickle
with open('protein_dict.pickle', 'rb') as f:
    protein_dict = pickle.load(f)

def PROT(string=None, filename=None, save=None):
    if string is not None:
        pass
    if filename is not None:
        with open(filename, 'r') as file:
            string = file.read()
    
    if len(string)%3 != 0:
        raise ValueError('dna length is not invalid')
    
    protein_string = ''
    for i in range(0, len(string), 3):
        if protein_dict[string[i:i+3]] == 'Stop':
            pass
        else:
            protein_string += protein_dict[string[i:i+3]]

    if save is not None:
        with open(save, 'w') as file:
            file.write(protein_string)

    return protein_string


def SUBS(string=None, substring=None, filename=None, save=None):
    if string is not None and substring is not None:
        pass
    elif filename is not None:
        with open(filename) as f:
            string, substring = f.read().splitlines()
    else:
        raise ValueError('input is invalid')

    k = len(substring)
    indexes = []
    for i in range(0, len(string)-k):
        if string[i:i+k] == substring:
            indexes.append(i+1)
    
    if save is not None:
        line = ''
        for i in indexes:
            line += ' '
            line += str(i)
        with open(save, 'w') as file:
            file.write(line[1:])

    return indexes


import numpy as np

def CONS(filename=None, save=None):
    sequences, headers = Fasta_read(filename)

    matrix = np.array([list(string) for string in sequences]).T
    nucleo = np.array(['A', 'C', 'G', 'T'])
    profile = np.array([np.sum(matrix == n, axis=1) for n in nucleo])

    c = nucleo[profile.argmax(axis=0)]
    c = ''.join(c)

    if save is not None:
        with open(save, 'w') as f:
            f.write(c)
            f.write('\n')
            f.write('A: ')
            f.write(' '.join(str(i) for i in profile[0]))
            f.write('\n')
            f.write('C: ')
            f.write(' '.join(str(i) for i in profile[1]))
            f.write('\n')
            f.write('G: ')
            f.write(' '.join(str(i) for i in profile[2]))
            f.write('\n')
            f.write('T: ')
            f.write(' '.join(str(i) for i in profile[3]))

    return c, profile



# RNA(filename='rosalind_rna.txt', save='Example_outputs/rosalind_rna_output.txt')
# REVC(filename='Example_datasets/rosalind_revc.txt', save='Example_outputs/rosalind_revc_output.txt')
# DNA(filename='Example_datasets/rosalind_dna.txt', save='Example_outputs/rosalind_dna_output.txt')
# print(FIB(31, 4))
# print(GC('Example_datasets/rosalind_gc.txt'))
# print(HAMM(filename='Example_datasets/rosalind_hamm.txt'))
# print(IPRB(30, 27, 28))
# print(PROT(filename='Example_datasets/rosalind_prot.txt', save='Example_outputs/rosalind_prot_output.txt'))
# print(SUBS(filename='Example_datasets/rosalind_subs.txt', save='Example_outputs/rosalind_subs_output.txt'))
# print(CONS(filename='Example_datasets/rosalind_cons.txt', save='Example_outputs/rosalind_cons_output.txt'))
