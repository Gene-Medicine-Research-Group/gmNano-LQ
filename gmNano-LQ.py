
#!/usr/bin/python

# gmNano-LQ.py
# v01 2025-03
# Steve Hyde

# The intended use of gmNano-LQ is to take a fastq file from a Nanopore plasmid DNA sequencing study and extract the sequence_ID, sequence length, Phred quality score
# This data can then be manipulated in eg Excel/Prism to display sequence length vrs quality which can be a way of spotting populations of different length eg IS element insertions

# gmNano-LQ takes as input a fastq file which is expetced to contain MANY DNA sequence reads with their associated Phred quality scores

# Nanopore services typiclly provide "assembly reads" as a .gz file to create the needed fastq use gunzip
# gunzip filename.fastq.gz
# nb this deletes the starting file, so do this on a copy!

# gmNano-LQ outputs either to the terminal or to a tab deliminated TXT file three columns: sequnce_ID, sequence length, Phred quality score

# suggested command line usage and arguments:
# python gmNano-DQ.py -h                <<<< provides brief help on the arguments gmNano-LQ expects
# -i filename.fastq                     <<<< REQUIRED the fastq input file 
# -o filename.txt                       <<<< OPTIONAL output filename, doesn't have to have .txt extension, but is recommended
# -of filename.fastq                    <<<< OPTIONAL outputs a fastq file with sequence_ID appended to include -gmNano-LQ-Dnnn-Qqqq- where nnn=sequence length, qqq=average Phred quality score
# -dn nnnn                              <<<< OPTIONAL minimum DNA sequence length nnnn (rounds floating point input down to nearest integer) bases to include in output
# -dm nnnn                              <<<< OPTIONAL maximum DNA sequence length nnnn (rounds floating point input down to nearest integer) bases to include in output
# -qn qqqq                              <<<< OPTIONAL minimum average Phred quality score qqqq (integer or float) to include in output
# -qm qqqq                              <<<< OPTIONAL maximum verage Phred quality score qqqq (integer or float) to include in output
# -s                                    <<<< OPTIONAL suppress the sequence ID in the output
# -v                                    <<<< OPTIONAL verbose mode - provides more information about the analysis to the terminal

# python gmNano-DQ.py -i 613435901_002_barcode27.fastq -o myFile.txt -dn 100 -dm 10000 -v

# useful information on fastq files and Phred quality scores:
# https://en.wikipedia.org/wiki/FASTQ_format
# https://en.wikipedia.org/wiki/Phred_quality_score

# average Phred quality score is calculated according to the approach used by Nanoplot (convert quality scores to error probability, average, convert back)
# https://gigabaseorgigabyte.wordpress.com/2017/06/12/oxford-nanopore-basecall-quality-scores/
# https://gigabaseorgigabyte.wordpress.com/2018/08/30/update-on-oxford-nanopore-basecall-quality-scores/


#needed for math.log10()
import math 

# setup parsed commands
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument('-i', '--input', dest='myInputFileName', type=str, help='Input file name to process')
parser.add_argument("-i", "--input", help="Input file name to process")
parser.add_argument("-o", "--output", help="Output analysis file name to write")
parser.add_argument("-of", "--output_fastq", help="Output fastq filename with sequence_id updated with sequence length and average Phred quality score")
parser.add_argument("-dn", "--DNAmin", help="Minimum DNA sequence length to include")
parser.add_argument("-dm", "--DNAmax", help="Maximum DNA sequence length to include")
parser.add_argument("-qn", "--Qualitymin", help="Minimum average Phred quality score to include")
parser.add_argument("-qm", "--Qualitymax", help="Maximum average Phred quality score to include")
#parser.add_argument("-s", "--sort", help="Sort the output file by sequence length", action="store_true")
parser.add_argument("-s", "--suppress_sequenceID", help="Don't provide sequence ID inforamtion in output", action="store_true")
parser.add_argument("-v", "--verbose", help="Provide greater explanation", action="store_true")
args = parser.parse_args()

if args.input:
    myInputFileName = args.input
    if args.verbose:
        print ('gmNano-LQ processing '+ myInputFileName)
else:
    print ("gmNano-LQ requires a fastq file. None provided. try --input filename.fastq")
    exit (1)

if args.output:
    myOutputFileName = args.output
    if args.verbose:
        print ("Will output results to "+myOutputFileName)
    myPrintFile = False
else:
    myPrintFile = True

if args.output_fastq:
    myOutputFASTQFileName = args.output_fastq
    if args.verbose:
        print ("Will output modified fastq file as "+myOutputFASTQFileName)
    myOutputFASTQ = True
else:
    myOutputFASTQ = False

if args.DNAmin == None:
    myMinSequenceLength = -1
else:
    myMinSequenceLength = math.floor(float(args.DNAmin))
    if args.verbose:
        print ("Minimum DNA sequence length to include "+str(myMinSequenceLength)+" bases")

if args.DNAmax == None:
    myMaxSequenceLength = 9999999999999
else:
    myMaxSequenceLength = math.floor(float(args.DNAmax))
    if args.verbose:
        print ("Maximum DNA sequence length to include "+str(myMaxSequenceLength)+" bases")

if myMaxSequenceLength < myMinSequenceLength:
    print ("gmNano-LQ requires a minimum DNA sequence length to be shorter or equal to the maximum DNA sequence length")
    exit (1) # exit with error

if args.Qualitymin == None:
    myMinQualityAverage = -1
else:
    myMinQualityAverage = float (args.Qualitymin)
    if args.verbose:
        print ("Minimum average Phred quality score to include "+args.Qualitymin)

if args.Qualitymax == None:
    myMaxQualityAverage = 1000
else:
    myMaxQualityAverage = float (args.Qualitymax)
    if args.verbose:
        print ("Maximum average Phred quality score to include "+args.Qualitymax)

if myMaxQualityAverage < myMinQualityAverage:
    print ("gmNano-LQ requires a minimum average Phred quality score to be less than or equal to the maximum Phred quality score")
    exit (1) # exit with error

if args.suppress_sequenceID:
    myOutputSequenceID = False
    if args.verbose:
        print ("Not including sequence ID in output")
else:
    myOutputSequenceID = True
    if args.verbose:
        print ("Including sequence ID in output")

if args.verbose and myPrintFile:
    print ("No output filename provided, will list output in terminal")

# open the input fastq file, count the number of sequence entries - each one spans 4 lines. Then close the file to get back to the top
myInputFile = open(myInputFileName, "r")
myTotalSequencesCount = int (sum(1 for _ in myInputFile) / 4)
myInputFile.close()

if args.verbose:
    print (str(myTotalSequencesCount)+" sequence entries found")

#open the input and output files for action
myInputFile = open(myInputFileName, "r")
if not myPrintFile:
    myOutputFile = open(myOutputFileName, "a")
if myOutputFASTQ:
    myOutputFASTQFile = open(myOutputFASTQFileName, "a")

while myTotalSequencesCount > 0:
    myTotalSequencesCount = myTotalSequencesCount -1

    myLine = (myInputFile.readline()) # read the sequence_id line
    mySequenceName = myLine.strip() #strip the newline from the sequence_id
    myLine = (myInputFile.readline()) #read the associated DNA base sequence
    mySequenceString = myLine.strip()
    mySequenceLength = str(len(myLine.strip())) # strip the final newline and determine length of DNA sequence 
    myLine = (myInputFile.readline()) # in fastq this line is always a "+" check/test - exit if not
    myCheckfastqFormat = myLine.strip()
    if not myCheckfastqFormat == "+":
        print ("Input file corrupt or not fastq format")
        exit (1)
    myLine = (myInputFile.readline()) # read the Phred quality score for each base
    myPhredString = myLine.strip() #strip the newline from the Phred quality string

    myRunningPhredProbability = 0 #variable to add up all the Phred score probabilities

    for character in myPhredString:
        myBasePhredScore = ord(character) - ord ('0') +15 #convert character to Phred score 0-50
        myBasePhredProbability = pow(10,-myBasePhredScore/10) #convert the Phred score to a probability
        myRunningPhredProbability = myRunningPhredProbability + myBasePhredProbability

    myAveragePhredProbability = myRunningPhredProbability/len(myPhredString)
    myAveragePhredScore = -10 * math.log10 (myAveragePhredProbability)

    mySequenceLengthAsInteger = int(mySequenceLength)


    if myPrintFile: #print the analysis at the terminal
        if ((mySequenceLengthAsInteger >= myMinSequenceLength) and (mySequenceLengthAsInteger <= myMaxSequenceLength) and (myAveragePhredScore >= myMinQualityAverage) and (myAveragePhredScore <=myMaxQualityAverage) ):
            if myOutputSequenceID:
                print (mySequenceName, sep='', end='\t')
            print (mySequenceLength, str(myAveragePhredScore), sep='\t')
    else: #output an analysis txt file
        if ((mySequenceLengthAsInteger >= myMinSequenceLength) and (mySequenceLengthAsInteger <= myMaxSequenceLength) and (myAveragePhredScore >= myMinQualityAverage) and (myAveragePhredScore <=myMaxQualityAverage) ):
            if myOutputSequenceID:
                myOutputFile.write(mySequenceName)
                myOutputFile.write("\t")
            myOutputFile.write(mySequenceLength)
            myOutputFile.write("\t")
            myOutputFile.write("{:.9f}".format(myAveragePhredScore))
            myOutputFile.write("\n")
            if args.verbose: # also output to terminal
                if myOutputSequenceID:
                    print (mySequenceName, sep='', end='\t')
                print (mySequenceLength, str(myAveragePhredScore), sep='\t')
    #end if myPrintFile:
    if myOutputFASTQ: #output a fastq file
        if ((mySequenceLengthAsInteger >= myMinSequenceLength) and (mySequenceLengthAsInteger <= myMaxSequenceLength) and (myAveragePhredScore >= myMinQualityAverage) and (myAveragePhredScore <=myMaxQualityAverage) ):
            mySequenceName = (mySequenceName+"-gmNano-LQ-D"+mySequenceLength+"-Q"+str("{:.9f}".format(myAveragePhredScore))+"-") #append -gmNano-LQ-Dnnn-Qqqq- where nnn=DNA sequence length qqq=average Phred quality score
            myOutputFASTQFile.write(mySequenceName) # not affected by -s as otherwise invalid fastq file format  
            myOutputFASTQFile.write("\n")
            myOutputFASTQFile.write(mySequenceString)
            myOutputFASTQFile.write("\n")
            myOutputFASTQFile.write("+")
            myOutputFASTQFile.write("\n")
            myOutputFASTQFile.write(myPhredString)
            myOutputFASTQFile.write("\n")
    #end if myOutputFASTQ:
#end while myTotalSequencesCount > 0:

myInputFile.close()

if not myPrintFile:
    myOutputFile.close()

if myOutputFASTQ:
    myOutputFASTQFile.close()

if args.verbose:
    print ("gmNano-LQ complete")

exit (0) # exit without error