#!/usr/bin/perl

# David Strubbe, October 2009
# See how many nodes of each type are available on the LRC supercluster
# by parsing the output of pbsnodes
# Usage: perl freenodes.pl

# list here the 'properties' to include
@prop_list = ("lawrencium", "debug-nano", "nano", "vulcan");

$pbsnodes_total = `pbsnodes`;
@pbsnodes = split('\n\n', $pbsnodes_total);

$bad_list = "";

foreach(@pbsnodes) {
   @node = split(' ', $_);
   $id = $node[0];
   $state = $node[3];
   $np = $node[6];
   $properties = $node[9];

   for ($iprop = 0; $iprop < @prop_list; $iprop++) {
       if($properties =~ $prop_list[$iprop]) {
           if($state =~ "down" || $state =~ "offline") {
               $bad_list .= `pbsnodes -ln $id\n`;
               next;
           }

           $total[$iprop][$np]++;
           $free[$iprop][$np]++;
           # this makes sure the element is created, even if none are free
           if(!($_ =~ "nusers=0,")) {
           $free[$iprop][$np]--;
           }
           # only entirely empty nodes are counted

           # without this, debug-nano nodes will be counted twice
           last;
       }
   }
}

format STDOUT_TOP =
Node Type     | Cores |  Empty | Total Up
==============|=======|========|=========
.

format STDOUT =
@<<<<<<<<<<<< | @>>>> | @>>>>> | @>>>>>>> 
$property $np $nfree $ntotal
.

for $iprop (0 .. @prop_list) {
   for $np (0 .. $#{$total[$iprop]}) {
       if($total[$iprop][$np] > 0) {
           $property = $prop_list[$iprop];
           $nfree = $free[$iprop][$np];
           $ntotal = $total[$iprop][$np];
           write;
       }
   }
}

print "\nNodes with problems:\n" . $bad_list;
