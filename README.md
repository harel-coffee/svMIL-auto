# integrative-omics
Variant annotation and prioritization

Running the code:

```
python main.py ../data/TPTNTestSet/TP.txt ../output/TP_annotations.txt
```

**IMPORTANT NOTE**

Files that are larger than 100 Mb (and ideally not much larger than 25 Mb) should not be uploaded to GitHub. This includes for example the output file, but also intermediate files. 

Also, the code does not run on the HPC using any of the default Python installations, it requires **Numpy 1.11.3** to work. I run the code using a Python virtual environment (virtualenv) with the versions shown below.

Versions on Marleen's laptop (working version):

**Python 2.7.13** <br />
**Numpy 1.11.3** <br />
**(optional) Neo4J 1.1.0b3 bolt driver** if you wish to switch the database type to Neo4J. The Neo4J feature is deprecated. <br />
