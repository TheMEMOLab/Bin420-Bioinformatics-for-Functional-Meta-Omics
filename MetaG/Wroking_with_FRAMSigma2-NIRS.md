# Wroking with FRAM Cluster Sigma2-NIRS

### What is this?

This document is intended to be a quick reference guide on the basic usage of the FRAM Sigma2 HPC cluster. For a complete reference please referer to the full documentation of [FRAM](https://documentation.sigma2.no/hpc_machines/fram.html).

## Login into FRAM

For login open a Command-line interfase (CLI) or Terminal  and type something like this. 

```bash
ssh auve@fram.sigma.no
```
> [!Important]
> Rember to change your user name

This will ask for your password, the one you got from Sigma2 by text. Type it

*Even you don't see anything the password is typed*

After the first attempt of login you will get something like:

```

Welcome to fram.sigma2.no

Documentation:      https://documentation.sigma2.no/
Support email:      support@nris.no
Request resources:  https://www.sigma2.no/apply-e-infrastructure-resources

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Latest news from:       https://opslog.sigma2.no

  o 2024-12-09: Saga maintenance stop
  o 2024-12-04: NIRD system software update
  o 2024-11-25: Reduced capacity due to NTNU cooling upgrade week 48
  o 2024-11-05: Slowness on Saga
  o 2024-10-30: Saga is partially down and very slow.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NOTE: The $USERWORK autocleanup period is at least 21 days and up to 42 days if
sufficient storage is available.

Project directories have backup as described in
https://documentation.sigma2.no/files_storage/backup.html

Last login: Sun Dec  8 18:07:26 2024 from ti0056q161-0632.bb.online.no
(miniforge3) [auve@login-2.FARM: ~]$
```
### FRAM main configuration 

**Named after the Norwegian arctic expedition ship Fram, the new Linux cluster hosted at UiT Arctic University of Norway is a shared resource for research computing capable of 1.1 PFLOP/s theoretical peak performance.

Fram is a distributed memory system which consists of 1004 dual socket and 2 quad socket nodes, interconnected with a high-bandwidth low-latency Infiniband network. The interconnect network is organized in an island topology, with 9216 cores in each island. Each standard compute node has two 16-core Intel Broadwell chips (2.1 GHz) and 64 GiB memory. In addition, 8 larger memory nodes with 512 GiB RAM and. The total number of compute cores is 32256.**

Let's take a look into this figure: 

![Cluster](https://github.com/TheMEMOLab/Bio326-NMBU/blob/main/images/cluster.png)

>[!Warning]
> **NEVER RUN A JOB IN THE LOGIN NODE!!! THE LOGIN NODE IS ONLY FOR LOOKING AND MANAGING FILES, INSTALL SOFTWARE AND WRITE SCRIPTS** 

**All users have access to the $HOME, so please DO NOT USE THE $HOME FOR STORAGE LARGE FILES (e.g. fastq, sam, databases). The $HOME directory is intended to allocate small software executables and SLURM scripts**

### Where can I storage large files? 

>[!Important]
> During the BIN420 Course all data, scripts, results and so must be written in
```
/cluster/projects/nn9987k
```

Let's go there and create a folder where you can work:


```bash
cd /cluster/projects/nn9987k
mkdir $USER
tree /cluster/projects/nn9987k/$USER
cd $USER
pwd
```

This last directory will be the directory we will use for all of the HPC sessions.

