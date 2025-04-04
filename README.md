# Homework 3: Parallelizing Genome Assembly
Due Date: April, 11th, 2025 at 11:59 PM.

This assignment is an introduction to writing distributed memory programs for applications with irregular communication patterns. Here we will implement a distributed hash table using UPC++. We will use our distributed hash table to evaluate one stage of a de novo genome assembly pipeline.

### Contributors: Brent, Kofi, Puja
## For Remote Students

Dear remote students, we are thrilled to be a part of your parallel computing learning experience and to share these resources with you! To avoid confusion, please note that the assignment instructions, deadlines, and other assignment details posted here were designed for the local students. You should check with your local instruction team about submission, deadlines, job-running details, etc. and utilize Moodle for questions. With that in mind, the problem statement, source code, and references should still help you get started (just beware of institution-specific instructions). Best of luck and we hope you enjoy the assignment!

## Background

Here is some background on de novo genome DNA assembly.  Strictly speaking, you don't need to understand all of this to complete the assignment, but it is interesting.  DNA assembly is the determination of the precise order of the nucleotides in a DNA molecule. A DNA molecule consists of four different bases, namely, adenine (A), guanine (G), cytosine (C), and thymine (T). For the purposes of this assignment we consider a DNA molecule to be a DNA strand, e.g. CTAGGAGCT (although in reality the number of bases for an actual DNA strand is much larger -- on the order of billions). De novo genome assembly is the assembly of strands of DNA in new genomes, without the use of a reference genome.

Unfortunately, we can not read the whole DNA strand in one go and therefore researchers have invented alternative methods to construct a genome. One such method is called shotgun sequencing. In this method, many copies of the original DNA strand are made.  Each copy is then fragmented randomly into pieces.  We cannot control how the copies are fragmented; they are randomly split into short contiguous fragments.  Next, we read each small fragment and each result is called a short read. Finally, we use all the short reads to reconstruct the original DNA strand. Figure 1 shows the process.

Figure 1: Shotgun Sequencing (source: http://people.mpi-inf.mpg.de/~sven/images/assembly.png).

The short reads returned from shotgun sequencing are significantly shorter than the actual DNA strand from which they were generated, so we must somehow combine them to create one long strand.  In addition, these short reads include sequencing errors, which makes assembling these fragments more difficult. There are methods to preprocess the short reads and remove the errors;  however, these methods are outside of the scope of this homework assignment and we refer the interested reader to [1].

The output of this first preprocessing stage of the de novo genome assembly pipeline is a set of unique DNA sequence fragments of length k, which we call k-mers. K-mers represent error-free DNA segments.  Each k-mer is associated with a forward and backward extension. A k-mer's backward extension is the base that precedes the k-mer in the DNA sequence, and its forward extension is the base that follows the k-mer in the DNA sequence.

Given a set of unique k-mers, we can build a special graph which is a compact representation of the connectivity among these k-mers. This special type of graph is called a de Bruijn graph, and in general a de Bruijn graph is used to represent overlaps between sequences of symbols.

In our particular case, the vertices of our de Bruijn graph are k-mers.  The edges of our de Bruijn graph represent k-mers which overlap by k-1 bases.  Since DNA is linear, each vertex in our de Bruijn graph is guaranteed to have at most two neighbors. Additionally, each vertex in the de Bruijn graph is unique since the k-mers are unique.  An example of such a de Bruijn graph is shown in Figure 2, where we illustrate a graph with k = 3. In this particular graph, nodes ATC and TCT are connected with an edge because they overlap in 2 bases (TC).

Figure 2: A de Bruijn graph with k = 3. We can identify three connected components (contigs) and six start nodes: GAT, TGA, AAT, TGC, AAC, CCG.  The contigs represented in this graph via connected components are: GATCTGA, AACCG and AATGC.

After building the de Bruijn graph, we can traverse it and find connected components called contigs.  Note that these connected components should be linear chains due to the nature of DNA.  Contigs are error-free DNA sequences significantly longer than the original reads. In later stages of the assembly pipeline, contigs will be linked together by leveraging information from the original reads to eventually produce a set of scaffolds which constitute the final result. In this assignment, we will be focusing on the parallel construction and traversal of the de Bruijn graph of k-mers (with a few assumptions and simplifications) and we refer the interested reader to [1] for more information on other stages of the pipeline.

## Problem Statement

The input to our problem is a set of unique k-mers, which are sequences of DNA bases of length k.  Each k-mer is associated with a forward extension, which is the base which precedes the k-mer in the DNA sequence, and a backward extension, which is the base that follows the k-mer in the DNA sequence.

Our k-mers are guaranteed to be unique and to overlap one another by exactly k-1 bases.  Our goal is to traverse these k-mers to produce a series of contigs, which are longer contiguous segments of DNA.  For a particular k-mer, we can determine the next k-mer in the sequence by taking the last k-1 bases and appending the forward extension.  In Python slicing notation, this is next_kmer = kmer[1:] + forward_extension.

Some k-mers are special, in that they have the special forward or backward extension F (which is not one of the regular A,T,C, or G bases we'd expect in DNA).  The special base F indicates the beginning or end of a read.  We call k-mers with the backward extension F start k-mers, since they are the k-mers which start a contig.  Similarly, k-mers with the forward extension F are end k-mers, since they end a contig.

### K-Mer Traversal Algorithm

To generate contigs, we will begin by inserting all the k-mers into a hash table.  While we are inserting the k-mers, we will check for any special start k-mers with backward extension F and append them to a list of start k-mers.  After building the hash table and identifying our start k-mers, we can generate a contig for each start k-mer by searching for the next k-mer in the sequence until we find an end k-mer with forward extension F, which signifies the end of the contig.

Here's what this looks like in Python pseudocode.

```
start_kmers = list()
kmer_hash = dict()

kmers = muh_kmer_dataset()

# Build hash table, collect start k-mers.
for kmer in kmers:
    # Suppose kmer.kmer is the k-mer string and
    # kmer is an object representing the k-mer,
    # including the forward and backward extension.
    kmer_hash[kmer.kmer] = kmer
    if kmer.backward_ext == 'F':
        start_kmers.append(kmer)

contigs = list()

# Traverse contigs, each starting with a start k-mer.
for start_kmer in start_kmers:
    contig = list()
    contig.append(start_kmer)
    while contig[-1].forward_ext != 'F':
        next_kmer = contig[-1].next_kmer()
        contig.append(kmer_hash[next_kmer])
    contigs.append(contig)
```

## Your Job

In the starter code, you will find a serial implementation of the contig generation stage of the de novo genome assembly pipeline. Your job is to parallelize this code using UPC++.

## Starter Code

The starter code is available on github at https://github.com/Berkeley-CS267/hw3  

Here's a rundown of the files you'll find in the starter code tarball.  kmer_hash.cpp and hash_map.hpp are the main files you'll be modifying.  You won't necessarily need to understand how the other source files work, just how to use the kmer_pair and pkmer_t data structures.

```
hw3 (root folder)
  |----CMakeLists.txt   - Builds the kmer_hash binary, which does contig generation.
  |----kmer_hash.cpp    - The main kmer hash file. /* YOU CAN MODIFY THIS */
  |----hash_map.hpp     - The hash table.  This is the main data structure you'll be parallelizing. /* YOU CAN MODIFY THIS */
  |----kmer_t.hpp       - Data structure which defines a kmer_pair, which is a k-mer and its extensions.
  |----pkmer_t.hpp      - Data structure which contains a k-mer as a packed binary string.
  |----packing.hpp      - Code that packs a k-mer character string into a binary string.
  |----read_kmers.hpp   - Code which reads k-mers from a file.
```

## Implementing a Hash Table

In hash_map.cpp, you'll find a serial implementation of a hash table using open addressing.  If this sounds unfamiliar to you, go ahead and skim the Wikipedia pages on hash tables and open addressing.  A hash table uses a hash function to decide at which location in an array to store an item.  This allows fast random access using a key value.  When two keys produce the same hash value, there is a conflict, since we cannot store two items at the same location in an array.  One way to resolve these collisions is open addressing, which resolves conflicts by probing the array with some deterministic probing strategy until a free spot is found.  The starter code uses open addressing with linear probing, which means we just try the next slot in the array until we find a free slot for our item.

### In Distributed Memory

In order to parallelize the code, you'll need to modify the data and used data members of the hash map to refer to distributed objects and then modify the insert() and find() methods accordingly. A common way to implement a distributed array in UPC++ is to create a C++ vector of upcxx::global_ptrs that point to arrays allocated in the shared segment on each rank.  You can then view the distributed array as one logically contiguous array and use arithmetic to write to the appropriate location in shared memory.

### [IMPORTANT] 

It is possible to write a trivial distributed hash table in UPC++ by using remote procedure calls to insert() or find() an item in a local std::unordered_map or std::map.  Do not submit this as your parallel solution, or you will lose credit.

## Scaling Experiments on Perlmutter Nodes

### [IMPORTANT] Scaling Experiments

Once your parallel code is working and optimized, run scaling experiments to see how good your performance is. Use the test.txt and human-chr14-synthetic.txt datasets.

    Perform multinode experiments using 1, 2, 4, and 8 nodes with 60 tasks per node. Graph the runtime of your parallel program for the different datasets as you increase the number of nodes. Also, graph your strong scaling efficiency. Are there any strange jumps or trends in your graphs? Did some datasets perform better than others? Include these graphs and analyses in your report. 

        How the performance change when moving from 60 tasks per node to 64 tasks per node?

    Run intra-node experiments using 1 node and varying the number of ranks per node (don't just use powers of two tasks per node). Also, see the note right below about memory segment size.

### [IMPORTANT] UPC++ Shared Memory Segment Size

The staff will test your solution with the default GASNET_MAX_SEGSIZE and UPCXX_SEGMENT_MB setting, so your solution must work in this setting for nodes from 1 to 8.

In your report, you must also include intra-node scaling on 1 node and [1-64] tasks per node (You can go up to 128 if you want). In this case, note that if you run larger data sets on a single node, you may need to increase the size of the UPC++ shared memory segment. This will allow you to allocate more shared arrays and objects. You can do this by setting the environment variable UPCXX_SEGMENT_ MB to a larger value. You may also need to increase GASNET_MAX_SEGSIZE. Common values are GASNET_MAX_SEGSIZE=4G and UPCXX_SEGMENT_MB =256, but you can choose any setting for the intra-node scaling as we will not be testing it. 

The staff will only test inter-node scaling using the default memory segment sizes. However, students must perform intra-node experiments to include them in the report, and in this case, students may vary the segment size as needed (and indicate this in their report). Please make sure the inter-node solution works with 64 and 60 tasks per node for the default segment sizes.

## Building and Running the Code

Note that in order to compile and run the code, you'll first need to load the appropriate UPC++ modules.
```
source modules.sh
```

You should then be able to build and run the starter code.

```
[demmel@perlmutter cs267_hw3_2025]$ mkdir build
[demmel@perlmutter cs267_hw3_2025]$ cd build
[demmel@perlmutter build]$ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=CC ..
[demmel@perlmutter build]$ cmake --build .
...
(For 19-mers)[demmel@perlmutter build]$ salloc -N 1 -A mp309 -t 10:00 --qos=interactive -C cpu srun -N 1 -n 1 ./kmer_hash_19 [my_dataset]
(For 51-mers)[demmel@perlmutter build]$ salloc -N 1 -A mp309 -t 10:00 --qos=interactive -C cpu srun -N 1 -n 1 ./kmer_hash_51 [my_dataset]
```


[my_dataset] is the dataset of k-mers you'd like to assemble contigs for.  You can find k-mer datasets in /global/cfs/cdirs/mp309/cs267-spr2025/hw3-datasets/ or $CFS/mp309/cs267-spr2025/hw3-datasets/. Both refer to the same location.

Here's a rundown of the various datasets.
```
  |----human-chr14-synthetic.txt  - 51-mers, A large dataset, based on the human chromosome 14.
  |----test.txt                   - 19-mers, A smaller dataset for testing.
  |----smaller (a folder, smaller datasets only for debugging)
       |----small.txt             - 19-mers, A smaller dataset for testing (1000 contigs)
       |----little.txt            - 19-mers, A little smaller dataset for testing (100 contigs)
       |----verysmall.txt         - 19-mers, A very small dataset for testing (10 contigs)
       |----tiny.txt              - 19-mers, A tiny dataset for testing (1 contig)
```

### [IMPORTANT] Recompiling for Different Sizes of K-Mers

Note that to run on the larger dataset, human-chr14-synthetic.txt, you'll need to recompile your code to handle 51-mers instead of 19-mers.  Simply modify the KMER_LEN macro at the top of packing.hpp, then recompile.  You'll need to do the same thing in reverse if you later want to switch back to 19-mers.  If you try to run your code on a dataset with the wrong kind of k-mer, it will automatically exit with an error telling you to modify packing.hpp and recompile. Note: The new CMake file should already be able to generate two compiled files, one with the length set to 19 and one with the length set to 51. You shouldn't have to manually set it anymore.

## Optimizing File I/O

File I/O between the project directory (where we've stashed the datasets) and the compute nodes can be quite slow, particularly for large files.  We recommend you (1) copy the datasets to a folder in your scratch space and (2) set the folder with your datasets to be striped.  Striping will spread out the files in your directory so that different parts are located on different physical hard disks (in Lustre file system lingo, these are called "OSTs", or object storage targets).  This can sometimes slow down serial I/O slightly, but will significantly increase I/O performance when reading one file from multiple nodes.
```
[demmel@perlmutter:$SCRATCH/cs267_hw3_2023]$ mkdir my_datasets
[demmel@perlmutter:$SCRATCH/cs267_hw3_2023]$ lfs setstripe -c 72 -S 8M my_datasets/
[demmel@perlmutter:$SCRATCH/cs267_hw3_2023]$ cd my_datasets
[demmel@perlmutter:$SCRATCH/cs267_hw3_2023/my_datasets]$ cp -r /global/cfs/cdirs/mp309/cs267-spr2025/hw3-datasets/* .
```

You can then call kmer_hash with the datasets in $SCRATCH/cs267_hw3_2023/my_datasets, and your runs should be somewhat faster.  Note that this is optional, and will not improve your timed performance, since we don't time file I/O in kmer_hash.  However, your runs will finish faster.  If you like, you can read more about file system performance here.

## Testing Correctness

You'll need to test that your parallel code is correct.  To do this, run your parallel code with the optional test parameter.  This will cause each process to print out its generated contigs to a file test_[rank].dat where [rank] is the process's rank number.  To compare, just combine and sort the output files, then compare the result to the reference solutions located in the same directories as the input files.
```
[demmel@perlmutter build]$ salloc -N 1 -A mp309 -t 10:00 -q debug --qos=interactive -C cpu srun -N 1 -n 32 ./kmer_hash_19 my_datasets/test.txt test
[demmel@perlmutter build]$ cat test*.dat | sort > my_solution.txt
[demmel@perlmutter build]$ diff my_solution.txt my_datasets/test_solution.txt
```

If diff prints a bunch of output telling you differences between the files, there's an issue with your code.  If it's quiet, your code is correct.  You can also use tools like md5sum or shasum to check whether you solution is correct.  Note that you should remove your output files (rm test*.dat) between test runs.

## Submission Details

Supposing your custom group name is XYZ, follow these steps to create an appropriate submission archive:

    Ensure that your write-up is located in your source directory. It should be named cs267XYZ_hw3.pdf

    [IMPORTANT] You must submit only your 'kmer_hash.cpp' and 'hash_map.hpp' files. ANY OTHER FILES SHOULD BE UNMODIFIED. If the code doesn't compile, you'll receive an email asking to resubmit a compiling/working version.

    From your build directory, run:
```
student@login005:~/hw3/build> cmake -DGROUP_NAME=XYZ ..
student@login005:~/hw3/build> cmake --build . --target package
```

This second command will fail if the PDF is not present.

    Confirm that it worked using the following command. You should see output like:

```
student@login005:~/hw3/build> tar tfz cs267XYZ_hw3.tar.gz 
cs267XYZ_hw3/cs267XYZ_hw3.pdf
cs267XYZ_hw3/kmer_hash.cpp 
cs267XYZ_hw3/hash_map.hpp 
```

    Submit your .tar.gz through bCourses.

## Write-up Details

Your write-up should contain:

        The names of the people in your group and each member's contribution.

        Graphs and discussion of the scaling experiments (inter-node using 64 tasks per node and intra-node varying number of tasks per node).

        Describe your implementation--how is your hash table organized?  What data structures do you use?

        Any optimizations you tried and how they impacted performance.

        A discussion of how using UPC++ implemented your design choices. How might you have implemented this if you were using MPI? If you were using OpenMP?

        [Optional] How does your solution perform compared to a straightforward solution using remote procedure calls to insert() or find() an item in a local std::unordered_map or std::map?

## Resources

### UPC++

    UPC++ training site: https://upcxx.lbl.gov/training.

    Guide HTML: https://upcxx.lbl.gov/docs/html/guide.html.

    Guide PDF: https://upcxx.lbl.gov/wiki/docs/guide.pdf.

    Spec PDF: https://upcxx.lbl.gov/wiki/docs/spec.pdf.

### Genome Assembly

    For general information about Evangelos Georganas' work on parallel genome assembly, upon which this homework is based, the first few chapters of his thesis are excellent.

    For information about the production version of this phase of the genome assembly pipeline, see this SC14 paper.

### Parallel Programming Languages

    You can read a high-level overview of UPC++ in this extended abstract.

    Titanium, UPC, and UPC++ are all parallel languages developed here in Berkeley.

    I find this early paper on Titanium quite interesting, as there are a lot of pioneering features for such an early project.  Smart pointers, in distributed memory, in 1998!

## Troubleshooting

    If you are having RPC calls randomly hanging/getting lost, try running export GASNET_OFI_RECEIVE_BUFF_SIZE=single.


    If you are still running into issues, per the FAQ on berkeleylab / upcxx / wiki / FAQ — Bitbucket 

> "Can you make an RPC call from within a RPC callback? When I try to wait on the nested RPC I get an error/hang."

Out of all the implementations that hang, this is the most likely reason. The simplest way is to avoid nested rpc calls. 

A simple approach that doesn't hang is to simply return the first rpc, and check the success status. If it fails, then invoke the next rpc from the master process. What this means is that suppose process 0 made the insert request to process 1, and process 1 no longer has any spots left. Instead doing a nested rpc call from process 1 to some process, say process 2, simply report failure back to process 0, and let process 0 make the rpc call to process 2.

## References

[1] Jarrod A. Chapman, Isaac Ho, Sirisha Sunkara, Shujun Luo, Gary P. Schroth, and Daniel S. Rokhsar. Meraculous: De novo genome assembly with short paired-end reads. PLoS ONE, 6(8):e23501, 08 2011.

## Additional References

https://people.eecs.berkeley.edu/~aydin/sc15_genome.pdf

https://upcxx.lbl.gov/docs/html/guide.html#remote-procedure-calls

https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf

Why are de Bruijn graphs useful for genome assembly?
https://pmc.ncbi.nlm.nih.gov/articles/PMC5531759/

Parallelizing Irregular Applications for Distributed Memory Scalability: Case Studies from
Genomics
https://escholarship.org/content/qt1400c4rh/qt1400c4rh_noSplash_a3cca84fb0d7680a6f6809fecedc182f.pdf

https://bitbucket.org/berkeleylab/upcxx/src/master/example/prog-guide/dmap.hpp

UPC++ v1.0 Training Materials: THIS!!!

https://bitbucket.org/berkeleylab/upcxx/wiki/Training

https://bitbucket.org/berkeleylab/cuf23/downloads/cuf23-upcxx.pdf

ExaBiome:
https://sites.google.com/lbl.gov/exabiome/home?authuser=0

https://bitbucket.org/berkeleylab/mhm2/src/master/

https://bitbucket.org/berkeleylab/upcxx/wiki/docs/system/perlmutter