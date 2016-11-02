#!/bin/bash
bzr export /tmp/ASVMLearning.zip
scp /tmp/ASVMLearning.zip webadmin@lasa.epfl.ch:/var/www/vhosts/asvm.epfl.ch/docs/downloads
