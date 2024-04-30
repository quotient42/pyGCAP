import subprocess

TAXON = input("<< Enter TAXON value: ")
subprocess.run(['bash', './src/init.sh', TAXON])