#!/bin/bash

#tentakel /bin/rm -f /tmp/*.*
tentakel /bin/rm -rf /tmp/bgovi/
#pdsh -g amino script
#tentakel /bin/rm -rf /tmp/bgovi/pdb

#tentakel /bin/mkdir -p /tmp/bgovi

#for x in `seq  0 7`
#do
#    for y in `seq  0 7`
#    do
#        node="compute-$x-$y"
#        echo $node
#        /usr/bin/ssh $node "rm -f /tmp/guerler/* < /dev/null >& /dev/null &"
#        /usr/bin/ssh $node "rm -f /tmp/* < /dev/null >& /dev/null &"
#    done
#done
