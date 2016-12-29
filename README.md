# wipo-st25-parse README File

This is a module for Biopython to allow for parsing and importing WIPO ST.25 sequence files as SeqRecord objects

# Objective

I had the misfortune of having to do some analysis on a bunch of patent sequence files in this WIPO ST.25 format.

The software provided by the USPTO for this purpose (https://www.uspto.gov/patents-getting-started/patent-basics/types-patent-applications/utility-patent/checker/patentin) has been "enhanced" by .NET. It is not exactly amenable to fitting into a sequence analysis pipeline. It did not open ST.25 files without import error messages that caused it to abort. I did not attempt to fix it.

This is a hacked together parser for ST.25 sequence files to get them into biopython objects. It works for the most part but it's pretty ugly, and the format specification is vague on a bunch of important points so lots of these files have subtle differences. Your milage may vary. Use at your own risk. I don't think enough people will have to work with these types of file to justify pulling it back into biopython-master as a module in SeqIO.parse, and if you do work with these files I would suggest keeping a bottle of scotch in your desk drawer. 

I am not ambitious enough to script something that would actually write files in this format. It should be easier than writing a parser, but I just do not want to do it. 

#Requirements

-	Python 2.X, 3.X (This should be agnostic to python version but I only tested it on the most recent releases of 2 & 3)

-	BioPython 1.67

	https://github.com/biopython/biopython