# Data QC & Pre-processing Handbook

Summary

During this practical you will take an Illumina paired end dataset and complete the following quality control workflow:

1.  Assess the quality of paired end data using **FastQC**
2.  Use **Trimmomatic** to trim/remove poor quality bases/reads
3.  Use **Trimmomatic** to remove 3’ adapter sequences from reads
4.  Re-assess the data using **FastQC**

This practical was created by Matthew Blades (BBASH, University of Leicester, UK) who also wrote the original version of this handbook.


## Getting Started with this handbook
Each stage of the workflow will be described in terms of what is being done and the command required to run the software.  Commands will be in the following format:
```
>$ command you need to type to run the software
```

Where:

``` >$ ```  is the command prompt.

The actual command prompt on your terminal will be in a format something like this:
```
[corso@benode01 'your_folder']$
```
e.g in my case the command prompt looks like this:
```
[corso@benode01 chiara]$
```


However, to save space in this manual I have replaced the actual command prompt with this:
```
>$
```

The text in   ```  this gray boxes ```  is the actual text you need to type into the command window to run the software.


**NOTE:  After typing in a command hit return to execute that command**


A command is generally the name of the software to be used followed by a list of parameters and options for how the software should be run. **For example**, the command below will run the software cutadapt. The text after cutadapt i.e.  -a, -i and -o are all cutadapt parameters followed by their input values

```
>$ cutadapt -a ATGAATCTA -i data_in.fastq -o data_out.fastq
```

**NOTE: Please be aware of spaces between text in commands as they are important!**


At certain points during this practical questions will be asked and a blank space left for you to fill in your answers. Questions are written in green text.
## Step1: Copy and check the data files

The data we will be using today comes from a [publication](http://goo.gl/dcKbB) where they used next generation sequencing to investigate population genetics of *Vibrio cholera*, following a cholera outbreak in the aftermath of the Haitian earthquake in 2010. The accession number for the complete data set is SRA039806, with the paired end data files having the accession number SRR308665.

Paired end 1 data file         = paired_end1.fastq

Paired end 2 data file         = paired_end2.fastq

**NOTE: remember that raw paired end data is contained in two separate files, one file for each end.**

First of all, open a connection to HPC Bender.

```
>$ ssh -X corso@bender.igb.cnr.it
```

The first task is to copy the data in your directory so you can use them. Use the ‘cd’ (change directory) command to move to your directory

```
>$ cd your_directory
```

The command ``` pwd ``` which stands for ‘print working directory’ will tell you the path of your current directory.

```
>$ pwd
```

In my case the output from the ``` pwd ``` command is:

```
/home/corso/chiara/
```
In order to keep things tidy in your directory it is good practise to create separate folders for the different exercises. We will now create a **``` Data_QC ```** folder in scratch directory then use the command ``` cp ``` (copy) to copy the files into it.

The ``` mkdir ``` command, short for ‘make directory’ will create a directory with the name you specify.

```
>$ mkdir Data_QC
```

Then ``` cd  ```  to the **``` Data_QC ```** directory you just created.
```
>$ cd Data_QC
```
We will now use the ``` cp ``` command to copy the data files from their current location
``` /home/corso/varcall2016/day2 ```

to your location (the Data_QC directory you just created).  Once you have run the ``` cp```  command use the ```  ls ```  (list) command to check the files have been copied to your directory.  
You will notice that there is **a space then a full stop ``` . ``` ** at the end of the 3 ``` cp ``` commands below, this is vital as it tells the  ``` cp `` command that you would like the files copied to your current directory.  

If you replaced the  ``` . ```  with a path e.g. home/MyDirectory/ then the files would be copied there instead. Note that the **TruSeq3-PE-2.fa** is a file of adapter sequences that ****Trimmomatic**** will need for the adapter removal stage later.

```
>$ cp /home/corso/varcall2016/day2/paired_end1.fastq .
>$ cp /home/corso/varcall2016/day2/paired_end2.fastq .
>$ cp /home/corso/varcall2016/day2/TruSeq3-PE-2.fa .
>$ ls
```
---
>
> You can also copy all the files with a single command, do you know how to do this?
>
---

Fastq files are usually too large to open with text editors, but there are Unix commands to help you look at your data.  Use the ``` more ``` command to view your fastq files (note the fastq format as described in the presentation earlier).

```
>$ more paired_end1.fastq
```
To exit the ``` more ``` command and return to the command prompt press ``` q ```

Try the ``` head ``` command which displays just the first 10 lines of a file

```
>$ head paired_end1.fastq
```
---
>
> **Question** – There is also a Unix command called ‘tail’, given what the ‘head’ command does what do you think ‘tail’ will do? Try it on one of your data files.
>
---
Now that you have checked the fastq files are there and viewed them you can proceed with the workflow.



## Step 2: Assess the quality of the data using **FastQC**

**FastQC** is a quality control tool for high throughput sequence data.

Launch **FastQC** as follows:

```
>$ /opt/bio/FastQC/fastqc &
```

If it is too slow, exit FastQC on bender and use it locally. Remember to transfer the fastq files to your laptop too!

From the FastQC window click:

```
> File > Open
```

Then browse to the location of the fastq files (this should be the ‘Data_QC’ directory you just created in your directory, or where you have downloaded them on your laptop) and select both of the paired end fastq files then click ‘OK’. It might take a minute or two for **FastQC** to open the files and analyse them.

Work your way through the analysis modules on the left hand side of the **FastQC** window, using the **FastQC** [documentation](http://goo.gl/8nSMkq) familiarize yourself with what each module is showing. Pay particular attention to the modules below as they will help direct the downstream trimming and adapter removal steps
- Basic Statistics
- Per Base Sequence Quality
- Per Sequence Quality Scores
- Overrepresented Sequences


**NOTE : There is a characteristic drop in quality scores of the 3’ bases of the reads (we will trim these later). The overrepresented sequences will also be removed.**

---
>
>Question – what is the total number of sequences in each of the paired end fastq files? (hint – use the basic statistics module)
>
>paired_end1.fastq
>
>paired_end2.fastq
>
>Question – what type of encoding is used in the fastq files?
>
>
>Question – what is the length of the sequences in the fastq files?
>
-----

**Do not close the FastQC window as we will compare the original files to the ones we will produce after adapter removal, and quality filtering.**



## Step 3: Use **Trimmomatic** to remove adapter sequences

**Trimmomatic** is a java tool for performing a range of trimming tasks on Illumina paired end and single end read data. The documentation can be found at this [link](http://www.usadellab.org/cms/uploads/supplementary/**Trimmomatic**/**Trimmomatic**Manual_V0.32.pdf.)


Go back to the overrepresented sequences module of **FastQC**. This is where **FastQC** would tell you if a significant proportion (>1%) of your reads are contaminated with adapter sequences. As you can see from the ‘Possible Source’ column, **FastQC** has found a number of reads contaminated with TruSeq adapter sequences, so we will run **Trimmomatic** to remove the adapter sequences.

The command below will run **Trimmomatic** and remove any adapter sequence from the reads.


**NOTE: The Trimmomatic command is long and won’t fit on one line below, therefore it is split over multiple lines, you will notice that each line has a backslash (\) at the end of it. This is standard command line procedure to indicate that a long command is split over multiple lines and continues on the next line below. It also helps to simplify long commands, make them easier to read and find errors in.
Note that the fourth and last line DOES NOT have a backslash at the end, this tells the terminal that this is the end of the command and to run that command when you next hit return. Where backslashes are used make sure the backslash is the last character on the line, i.e. no spaces after it. Obviously, it is not mandatory to use backslashes and divide the command in several lines, but we are doing it here for clarity.**

```
$ java -jar /opt/bio/**Trimmomatic**.jar PE -phred33 -trimlog logfile \

paired_end1.fastq paired_end2.fastq Left_paired.fastq \

Left_unpaired.fastq Right_paired.fastq Right_unpaired.fastq \

ILLUMINACLIP:TruSeq3-PE-2.fa:2:40:15 MINLEN:36

```

The parameters used for **Trimmomatic** are defined as follows:
	PE	(data is paired end)

-phred33                             (Quality scores are 33 offset)

-trimlog logfile                    (name of logfile for summary information)

paired_end1.fastq              (name of input fastq file for left reads)

paired_end2.fastq              (name of input fastq file for right reads)

Left_paired.fastq                (paired trimmed output fastq file for left reads)

Left_unpaired.fastq            (unpaired trimmed output fastq file for left reads)

Right_paired.fastq              (paired trimmed output fastq file for right reads)

Right_unpaired.fastq          (unpaired trimmed output fastq file for right reads)

ILLUMINACLIP                  (parameters for the adapter clipping)

TruSeq3-PE-2.fa     (text file of adapter sequences to search for)

:2:40:15                   (adapter-read alignment settings – see manual)

MINLEN:36                         (delete reads trimmed below length MINLEN)

---
>Question – According to the **Trimmomatic** screen output, what is the number and percentage of read pairs that ‘both survived’ adapter trimming?

>Number
>Percentage


>Given that the original fastq files each contained 4052587 reads, how many pairs of reads have been trimmed and then deleted by **Trimmomatic** in this step?

----

**NOTE:  Remember that Trimmomatic only deletes reads if the length after trimming of adapter sequences is less than MINLEN (which we set to 36bp).**







## Step 4: Use **Trimmomatic** to trim low quality bases

The **FastQC** ‘Per Base Sequence Quality’ and ‘Per Sequence Quality Scores’ modules have already told us that there could be some issues with the quality scores of the last few bases of the reads. We also know from the presentation that the quality of 3’ bases of sequence reads does tend to decrease. We will use **Trimmomatic** to trim poor quality bases from the 3’ end of the reads. **Trimmomatic** also checks the 5’end for poor quality bases. The command below will carry out the trimming on the adapter trimmed fastq files we created above.
```
$ java -jar /opt/bio/**Trimmomatic**.jar PE -phred33 -trimlog logfile2 Left_paired.fastq Right_paired.fastq Left_trim_paired.fastq Left_trim_unpaired.fastq Right_trim_paired.fastq Right_trim_unpaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

The parameters used for **Trimmomatic** are defined as follows:

PE	(data is paired end)

-phred33                            (Quality scores are 33 offset)

-trimlog logfile2                  (name of logfile for summary information)

Left_paired.fastq                (name of input adapter trimmed left fastq file)

Right_paired.fastq              (name of input adapter trimmed right fastq file)

Left_trim_paired.fastq        (paired trimmed output fastq file for left reads)

Left_unpaired.fastq            (unpaired trimmed output fastq file for left reads)

Right_paired.fastq              (paired trimmed output fastq file for right reads)

Right_unpaired.fastq          (unpaired trimmed output fastq file for right reads)

LEADING:3                        (Trim 5’ bases with quality score < 3)

TRAILING:3                       (Trim 3’ bases with quality score < 3)

SLIDINGWINDOW:4:15     (see manual for explanation)

MINLEN:36                        (delete reads trimmed below length MINLEN)


---
>
>Question – According to the **Trimmomatic** screen output, what is the number and percentage of read pairs that ‘both survived’ low quality base trimming?

>Number

>Percentage


>Given that the adapter trimmed fastq files contained 2165674 reads, how many pairs of reads have been trimmed and then deleted by **Trimmomatic** in this step?
>
---



## Step 5: Assess again the quality of the data using **FastQC**

Use **FastQC** to open the trimmed fastq files.
Go through the analysis modules on the left hand side to see what has changed compared to the original fastq files.

---
>
>Question – What are the main changes you see in the read quality?
>
---

**Well done……you have successfully QC’d and pre-processed a set of paired end read files!**
