---
course: NGS for evolutionary biologists: from basic scripting to variant calling
title: getting connected and organizing the space
author: Enza Colonna
credits: Mario Aversano 
time:
---

#Summary

- [Bender](#section-id-30)
  - [Get there](#section-id-32)
   - [Workspace](#section-id-99)
    - [GATEWAY](#section-id-103)
    - [WORKER NODES](#section-id-128)
    - [FILE SYSTEM](#section-id-156)
  - [Working on Bender](#section-id-177)
    - [1. Interactive mode](#section-id-179)
    - [2. Submitting jobs to a job scheduler](#section-id-197)
      - [a) Preparing the PBS script](#section-id-211)
      - [b) Submitting the job to PBS](#section-id-249)
      - [c) Checking the job status](#section-id-256)




# Bender 

<div id='section-id-32'/>

## Get there

We will be hosted for this course on a machine named Bender. This machine is on the top floor in this building and has all the software we need installed and tested.

To use remote machine the first step is to get connected to them and one protocol for connection is called [**Secure Shell**](https://en.wikipedia.org/wiki/Secure_Shell) or **SSH**. 
![bender](img/benderssh.png)
Before using the remote machine you should have talked to the machine administrator that will create an account for you with a username and a password. Once you obtain an account,  you can use SSH from your terminal. The basic instruction for connection is the command ```ssh``` followed by the username and the [IP address](https://en.wikipedia.org/wiki/IP_address) of the machine. Some times to  simplify we use a literal synonim of the IP address: 

```
$ ssh username@111.111.1111.11

$ ssh usrname@machinename

```

Once connected, if it is the first time you will be asked to confirm that you really want to connect to the machine, otherwise you will be just asked to type a password. 

*Beware: you won't se any letter on your screen while you type the password!*


We have created a shared user name for this course that is ```corso```, while the domani name for the machine is ```bender.igb.cnr.it```. Therefore to  **connect to Bender** we will type: 

```
auser@itslaptop:$ ssh -X corso@bender.igb.cnr.it
      __      __       .__
      /  \    /  \ ____ |  |   ____  ____   _____   ____
      \   \/\/   // __ \|  | _/ ___\/  _ \ /     \_/ __ \
       \        /\  ___/|  |_\  \__(  <_> )  Y Y  \  ___/
        \__/\  /  \___  >____/\___  >____/|__|_|  /\___  >
             \/       \/          \/            \/     \/
        __           __________                   .___
      _/  |_  ____   \______   \ ____   ____    __| _/___________
      \   __\/  _ \   |    |  _// __ \ /    \  / __ |/ __ \_  __ \
       |  | (  <_> )  |    |   \  ___/|   |  \/ /_/ \  ___/|  | \/
       |__|  \____/   |______  /\___  >___|  /\____ |\___  >__|
                             \/     \/     \/      \/    \/

[corso@bender ~]$ 


```
The ``` -X ``` option allow some graphical visualization. Note how username has canged at the prompt after ssh. 

<div id='section-id-45'/>



## Workspace

Bender is a very organized machine. Let's look a little bit closer to this:  

![bender](img/benderscheme.png)

<div id='section-id-103'/>

### GATEWAY 

This is the access point to Bender form the public network when we ssh to Bender we arrive here as first instance. *We are not allowed to work in this space*, it is instead necessary to go in one of the Worker Nodes (see below).  

<div id='section-id-128'/>

### WORKER NODES 

You can use ```ssh``` to connect to a worker node: 

```
MARIO 
```
but to make life easier we have created a shortcut and to go to node one you only need to type: 

```
[corso@bender ~]$ 01
Last login: Wed Apr 27 09:09:45 2016 from bender
[corso@benode01 ~]$ 

```


### FILE SYSTEM 

Once connected you will be sharing the same workspace with other users. It is therefore important to respect some rules and beware of other people folders! 

The file system is subdivided in these main folders: 

```
drwxrwxr-x  4 corso corso 4096 26 apr 12:34 chiara
drwxrwxr-x  2 corso corso 4096 27 apr 14:05 enza
drwxrwxr-x  3 corso corso 4096 22 apr 16:33 igv
drwxrwxr-x 23 corso corso 4096 27 apr 14:05 students
drwxrwxr-x  7 corso corso 4096 27 apr 14:00 varcall2016

```

- The ```varcall2016``` folder contains all the files that we will need for the practicals and it is subdivided in days and projects: 

```
[corso@benode01 ~]$ ls -l varcall2016/
totale 20
drwxrwxr-x 3 corso corso 4096 22 apr 17:19 day1
drwxrwxr-x 2 corso corso 4096 22 apr 15:48 day2
drwxrwxr-x 4 corso corso 4096 25 apr 10:40 day3
drwxrwxr-x 2 corso corso 4096 25 apr 12:44 day4
drwxr-xr-x 4 corso corso 4096 26 apr 14:03 project_1
drwxr-xr-x 4 corso corso 4096 26 apr 14:03 project_2

```

- Within ```students``` there is a folder with your name. **All your work should be done within this folder!** Every time check that you are in the right place. TO navigate to your folder: 

```
[corso@benode01 ~]$ cd students/myname

```




## Working on Bender 

<div id='section-id-179'/>

### 1. Interactive mode
Interactive mode is a command line shell which gives immediate feedback for each statement. The user can interact with the shell.

An example of  interactive use of the shell is:
```
echo thiscourseiscool
thiscourseiscool

```

A little more complicated example of interactive use of the shell is:
```
python myprogram.py inputfile1.txt inputfile2.someextension > myoutput.out

```
Working in interactive mode is OK for small tasks, but if the command line we are executing  takes a lot of memory we might get stuck, and what is worst we will use all the machine power preventing other people from using it.


<div id='section-id-197'/>

### 2. Submitting jobs to a job scheduler

Bender is a shared machine, that means that many users use it at the same time, therefore it is advised and polite to use queues for running commandlines.

![queue](img/queue.png)

Bender has a job scheduling system named PBS. To become familiar with PBS read its interesting story on the [PBS wiki page](https://en.wikipedia.org/wiki/Portable_Batch_System), go through the [official documentation](http://www.pbsworks.com/SupportGT.aspx?d=PBS-Professional,-Documentation), or google for one of the many tutorials available.

In the official documentation we read that PBS consists of a set of commands and system daemons/service that perform these tasks:

- **Queuing jobs**: PBS collects jobs (work or tasks) to be run on one or more computers. Users submit jobs to PBS, where they are queued up until PBS is ready to run them.
- **Scheduling jobs**: PBS selects which jobs to run, and when and where to run them, according to the policy specified by the site administrator.  PBS allows the administrator to prioritize jobs and allocate resources in a wide variety of ways, to maximize efficiency and/or throughput.
- **Monitoring jobs**: PBS tracks system resources, enforces usage policy, and reports usage.  PBS tracks job completion, ensuring that jobs run despite system outages.  

<div id='section-id-211'/>

#### a) Preparing the PBS script

To use PBS we need to  wrap all the relevant instruction in a text file. The PBS file is made of two parts:
- the PBS instructions
- the commandline

It looks like this:

```


#!/bin/bash
#PBS –N myjobname      #name of the job
#PBS -o job.out       #output file
#PBS -e job.err       #error file
#PBS -l select=1:ncpus=20:mpiprocs=20:mem=122GB  #resources
#PBS -l walltime=1:00:00                  #hh:mm:ss
#PBS -q <queue>              #chosen queue
#PBS -A <my_account>   #name of the account
#PBS -W group_list=<group>   #name of effective group for reservation

echo thiscourseiscool
```
Save the file with a reasonable name, e.g. `thisspecificjobpbs.sh`

In this file all the lines starting with `#PBS`  are options for the PBS scheduler, and `-N` in `#PBS -N` indicates the option for giving a name to the job. The other lines are specific to the task we are doing (e.g. printing "thiscourseiscool"). Finally, the file starts with a line `#!/bin/bash` that tells the machine that this is a shell script.

There are many PBS options and they all have default values. Sometimes we need to change the default values, for example we want to change the jobname at every job. Some other times  the default values are OK, and we don't need to include the option specification in the PBS file.

See here(http://www.hpc.cineca.it/sites/default/files/PBSProUserGuide13.0beta.pdf) for a list of PBS options.

A little clarification on the `-l` resources option:

- select = number of nodes requested
- ncpus = number of cpus per node requested
- mpiprocs = number of mpi tasks per node
- mem = RAM memory per node  


<div id='section-id-249'/>

#### b) Submitting the job to PBS

Jobs are submitted using:
```
qsub thisspecificjobpbs.sh
```

<div id='section-id-256'/>

####  c) Checking the job status

Once submitted, it is possible to check the job status using:
```
qstat
```

Most likely you will see many jobs running among which your(s).

```
vcolonna@bender.igb.cnr.it$ qstat
Job id            Name             User              Time Use S Queue
----------------  ---------------- ----------------  -------- - -----
51133.node001     evolfun          cdarwin                  0 Q parallel        
69083.node001     nicegenes        rfisher                  0 Q parallel        
77560.node001     beer             bender                   0 Q parallel        
165186.node001    nicepeas         gmendel           308:12:2 R parallel         
165866.node001    myjobname        vcolonna          00:00:00 E parallel  

```



