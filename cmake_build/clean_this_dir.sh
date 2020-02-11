
#!/bin/bash

shopt -s extglob

rm -vrf !("README"|"user_config"|"clean_this_dir.sh")
