# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
from math import sqrt
from Bio import SeqIO
from bioinfokit import analys


def countOccur(array, resultlist):
    file2 = open("result.txt", "a")
    for x in range(100):
        amount = 0
        for j in array:
            if x == int(j):
                amount = amount + 1

        resultlist.append(amount)
        file2.write(str(amount) + "\n")


def countAmount(sequence, frequency):
    length = len(sequence)
    cCount = sequence.count("C")
    gCount = sequence.count("G")

    Sum = cCount + gCount
    freq = round((Sum / length) * 100, 0)
    frequency.append(freq)

    file = open("counts.txt", "a")
    file.write(str(freq) + "\n")

    return freq


def calculate(sequencelist, frequency, resultlist):
    array = []
    for seq_record in SeqIO.parse("reads_for_analysis.fastq", "fastq"):
        sequence = seq_record.seq
        sequencelist.append(seq_record.id)
        sequencelist.append(sequence)
        array.append(countAmount(sequence, frequency))
    countOccur(array, resultlist)


def encodingCheck():
    analys.format.fq_qual_var(file="reads_for_analysis.fastq")


def peakSearch(frequency, index, answerlist, sequencelist):
    count = 0
    for value in frequency:
        if int(value) == index:
            if count < 5:
                answerlist.append(sequencelist[(frequency.index(value) * 2)])
                answerlist.append("\n")
                answerlist.append(sequencelist[(frequency.index(value) * 2)+1])
                answerlist.append("\n")
                sequencelist.pop(frequency.index(value) * 2)
                sequencelist.pop(frequency.index(value) * 2)
                frequency.remove(value)
                count = count + 1
    answerlist.append("\n")


sequences = []
frequencies = []
result = []
encodingCheck()
calculate(sequences, frequencies, result)

first = 34
second = 54
third = 70

answer = []

peakSearch(frequencies, first, answer, sequences)
peakSearch(frequencies, second, answer, sequences)
peakSearch(frequencies, third, answer, sequences)

file3 = open("sequences.txt", "a")

for x in answer:
    file3.write(str(x))
    print(x)
