from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt 
import sys 

# Program do wizualnego porównania par sekwencji DNA 

def readFastaFile(file):
    """
    Funkcja wczytująca plik FASTA
    Args:
        file (str): nazwa pliku FASTA
    Returns:
        sequence (str): zczytana sekwencja DNA
    """
    fastaSequence = SeqIO.parse(open(file),'fasta')
    for fasta in fastaSequence:
        sequence = str(fasta.seq)
    return sequence

def createMatrix(seq1, seq2):
    """
    Funkcja tworząca macierz na podstawie sekwencji DNA
    Args:
        dna1 (str): pierwsza sekwencja DNA
        dna2 (str): druga sekwencja DNA
    Returns:
        data: macierz
    """
    
    data = [[(seq1[i:i + 1] != seq2[j:j + 1]) \
        for i in range(len(seq1))] \
        for j in range(len(seq2))] 
    
    data = np.array(data).astype(int)
    data = np.where((data == 1) | (data == 0), data^1, data)
    return data
    
def createPlot(seq1, seq2, window, threshold, file1, file2):
    
    """
    Funkcja tworząca wykres kropkowy
        Args:
        dna1 (str): pierwsza sekwencja DNA
        dna2 (str): druga sekwencja DNA
        file1 (str): nazwa pliku/sekwencji
        file2 (str):nazwa pliku/sekwencji
        window (int): okno
        threshold (int): próg

    """
    
    data = createMatrix(seq1, seq2)
    suptitle = "Porównanie " + file1 + " i " + file2
    title = "Okno: " + str(window) + " i próg: " + str(threshold)
    ax = plt.subplot()
    
    if len(seq1) and len(seq2) < 12:
        plt.xticks(np.arange(len(list(seq1))),list(seq1))
        plt.yticks(np.arange(len(list(seq2))),list(seq2))
        ax.tick_params(top=True, bottom=False)
        ax.tick_params(labeltop=True, labelbottom=False)
        ax.xaxis.set_label_position('top')
    
    plt.imshow(data, cmap=plt.cm.Greys)
    plt.ylabel(file1)
    plt.xlabel(file1)                      
    plt.ylabel(file2)
    plt.title(title)    
    plt.suptitle(suptitle) 
    plt.savefig('dotplot.png')
    plt.show() 

    
    
def createFilterPlot(seq1, seq2, matrix, window, threshold, file1, file2):
    
    """
    Funkcja tworząca wykres kropkowy
        Args:
        dna1 (str): pierwsza sekwencja DNA
        dna2 (str): druga sekwencja DNA
        file1 (str): nazwa pliku/sekwencji
        file2 (str):nazwa pliku/sekwencji
        window (int): okno
        threshold (int): próg

    """
    
    suptitle = "Porównanie " + file1 + " i " + file2
    title = "Okno: " + str(window) + " i próg: " + str(threshold)
    ax = plt.subplot()
    
    if len(seq1) and len(seq2) < 12:
        plt.xticks(np.arange(len(list(seq1))),list(seq1))
        plt.yticks(np.arange(len(list(seq2))),list(seq2))
        ax.tick_params(top=True, bottom=False)
        ax.tick_params(labeltop=True, labelbottom=False)
        ax.xaxis.set_label_position('top')
    
    plt.imshow(matrix, cmap=plt.cm.Greys)
    plt.ylabel(file1)
    plt.xlabel(file1)                      
    plt.ylabel(file2)
    plt.title(title)    
    plt.suptitle(suptitle)
    plt.savefig('dotplot2.png')
    plt.show() 
    

def filterMatrix(matrix, window, threshold):
    
    """
    Funkcja filtrująca macierz
    Args:
        matrix: macierz
        window (int): okno
        threshold (int): próg
    Returns:
        matrix3: przefiltrowana macierz
    """
    
    matrix3 = matrix.copy()
    
    for i in range (-(len(matrix3)-1), len(matrix3[0])):
        diagonal = matrix3.diagonal(i)
        
        for j in range(0, len(diagonal)):

            if (i+window) > len(diagonal):
                if round(sum(diagonal[j:len(diagonal)-1])) < threshold:
                    if i < 0:
                        matrix3[-i+j,0+j] = 0
                    else:
                        matrix3[0+j,i+j] = 0
            else:
                if round(sum(diagonal[j:j+window])) < threshold:
                    if i < 0:
                        matrix3[-i+j,0+j] = 0
                    else:
                        matrix3[0+j,i+j] = 0
    return matrix3
    
    
# Wybór, czy chcemy wczytać plim FASTA, czy wprowadzić sekwencje DNA w konsoli
    
method = int(input("Wybierz nr metody:\n 1 - wprowadzenie danych w konsoli\n 2 - wczytanie pliku typu FATSA\n"))
stop = 0

# Wprowadzenie i przetworzenie danych dla metody 1

if method == 1:
    
    dna1 = str(input("Podaj pierwszą sekwencję DNA:\n"))
    for letter in dna1:
        if letter not in 'ATCG':
            print("Podano nieprawidłowe dane.")
            sys.exit()

        
    if stop == 0:
        dna2 = str(input("Podaj drugą sekwencję DNA:\n"))    
        for letter in dna2:
            if letter not in 'ATCG':
                print("Podano nieprawidłowe dane.")
                sys.exit()

    dnaName1 = "DNA 1"
    dnaName2 = "DNA 2"
    createMatrix(dna1, dna2)
    createPlot(dna1, dna2, 0, 0, dnaName1, dnaName2)


# Wprowadzenie i przetworzenie danych dla metody 2

elif method == 2:
    file1 = str(input("Podaj nazwę pierwszego pliku:\n"))
    file2 = str(input("Podaj nazwę drugiego pliku:\n"))
    dna1 = readFastaFile(file1)
    dna2 = readFastaFile(file2)
    createPlot(dna1, dna2, 1,  1, file1, file2)

else:
    print("Prosze wybrać prawidłowy numer.")
    

# Wybór dalszych czynności


nextStep = int(input("Wybierz co chcesz zrobić dalej:\n 1 - zmienić wartości okna i progu\n 2 - zakończ\n"))
 
if nextStep == 1:
    newWindow = int(input("Podaj wartość okna:\n"))
    
    if newWindow >= 1 and newWindow <= min(len(dna1), len(dna2)):
        newThreshold = int(input("Podaj wartość progu:\n"))
        
        if newThreshold >= 1 and newThreshold <= min(len(dna1), len(dna2)):
            matrix = createMatrix(dna1, dna2)
            filterMatrix(matrix, newWindow, newThreshold)
            createFilterPlot(dna1, dna2, filterMatrix(matrix, newWindow, newThreshold), newWindow, newThreshold, dnaName1, dnaName2)
            
        else:
            print("Podano nieprawidłowe dane.")
    
    else:
        print("Podano nieprawidłowe dane.")
        
    
elif nextStep == 2:
    sys.exit()
    
else:
    print("Niepoprawne dane.")