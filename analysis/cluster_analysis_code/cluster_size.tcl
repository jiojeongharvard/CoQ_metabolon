proc cluster_size {molid} {

package require pbctools

# get the file info of the top molecule
set name [molinfo $molid get filename]

# file name without tail 
set basename [file rootname [file tail [lindex $name 0]]]

# file path
set path [file dirname [lindex [lindex $name 0] 0]]

# get traj frame number
set numframes [molinfo $molid get numframes]

set outDir ${path}/cluster_results_skip_first_frame

set center_atom_name "index"


# UNIFORM, 50 each, COQ model
# Append all odd integers from 1 to 500 in system.data --> all even ints from 0 to 499 in VMD
for {set i 0} {$i <= 499} {incr i 2} {
    append center_atom_name " $i"
}

# Append all integers from 501 to 550/549 in system.data --> 500 to 549/548 in VMD
for {set i 500} {$i <= 549} {incr i} {
    append center_atom_name " $i"
}


set center_atom_sel [atomselect $molid $center_atom_name]
set center_atom_list [$center_atom_sel get index]

# if distance between two chains is larger than 2 rg, then they are not in the same cluster
set threshDist [expr (25)]

# dcd output frequency
set dcdFreq 2000

# timestep in simulation (fs)
set timestep 25

# set stride
set stride 1000

#######################################################
# Cluster progs
 
# Do two lists have any of the same elements?
proc checkOverlap {aList bList} {
    foreach a $aList {
    set i [lsearch -sorted -integer $bList $a]
    if {$i >= 0} {
        return 1
    }
    }
    return 0
}

proc searchOverlap {neighListList} {
    set n [llength $neighListList]
    set n1 [expr {$n-1}]
    # Loop over pairs of elements.
    for {set i 0} {$i < $n1} {incr i} {
    set li [lindex $neighListList $i]
    for {set j [expr {$i+1}]} {$j < $n} {incr j} {
        set lj [lindex $neighListList $j]
        if {[checkOverlap $li $lj]} {
        # If they overlap, merge the lists. 
        # If '-unique' is specified, then only the last set of duplicate elements found in the list will be retained.
        set lij [lsort -unique -integer [concat $li $lj]]
        # Delete the list $j.
        # Delete the later one since it's higher in the array, otherwise index of list $li will be affected.
        set neighListList [lreplace $neighListList $j $j]
        # Delete list $i and insert the merged list.
        set neighListList [lreplace $neighListList $i $i $lij]
        return $neighListList
        }
    }
    }
    return $neighListList
}

# Make unique clusters from lists of neighbors.
proc getClusters {neighListList} {
    set n [llength $neighListList]
    set n0 0
    while {$n != $n0} {
    set n0 $n
    set neighListList [searchOverlap $neighListList]
    set n [llength $neighListList]
    }
    return $neighListList
}

proc mean {l} {
    set sum 0.0
    foreach i $l {
    set sum [expr {$sum + $i}]
    }
    return [expr {$sum/[llength $l]}]
}


#######################################################
## main body


# Compute the interval in nanoseconds.
set frameInterval [expr {1e-6*$timestep*$dcdFreq}]

set displayPeriod 10

if {$stride < 1} { set stride 1 }

# Get the time change between frames in nanoseconds.
set dt [expr {$frameInterval*$stride}]

if 0 {
# Check if any of the files already exist
if {[file exists $outDir/clusterSizeMean.dat] || \
    [file exists $outDir/clusterNum.dat] || \
    [file exists $outDir/clusterDist.dat] || \
    [file exists $outDir/clusterSizeMax.dat] || \
    [file exists $outDir/clusterList.dat] || \
    [file exists $outDir/clusterMaxList.dat] || \
    [file exists $outDir/clusterCenterOfMass.dat] || \
    [file exists $outDir/clusterRadiusList.dat]} {
    error "One or more files already exist. Terminating script..."
}
}

# Open the output file.
set out [open $outDir/clusterSizeMean.dat w]
set outN [open $outDir/clusterNum.dat w]
set outDist [open $outDir/clusterDist.dat w]
set outMax [open $outDir/clusterSizeMax.dat w]
set outList [open $outDir/clusterList.dat w]
set outMaxList [open $outDir/clusterMaxList.dat w]

#set outComList [open $outDir/clusterCenterOfMass.dat w]
#set outRadList [open $outDir/clusterRadiusList.dat w]
#set outLigList [open $outDir/clusterLigList.dat w]


# Loop over the dcd files.
set nFrames0 0

# Move forward computing at every step, start with frame f = 1 not 0
for {set f 1} {$f < $numframes} {incr f} {
    # Get the time in nanoseconds for this frame.
    set t [expr {($nFrames0+$f)*$dt}]
    # Get the positions at the current frame.
    molinfo $molid set frame $f

    # Make a list of all of the neighbors of each residue.
    set neighListList {}

    foreach centeratom $center_atom_list {
    set s [atomselect $molid "($center_atom_name) and pbwithin $threshDist of index $centeratom"] 
# here i deleted 'and' compared with original code as this with make it faster -- see vmd short circuiting
    # here the 'and' means logic, rather than 'A + B'
    # test what's inside $s, selection of s is pretty important -- for residue x, s gives all center beads within the threshDist.
    # puts "[lindex $s]"
    
    # what if {1 3} {1 3 4} in neighListList? -- A: use getCluster function to merge two lists.
    set neighList [lsort -unique -integer [$s get residue]]
    lappend neighListList $neighList
    $s delete
    }

    # Get the clusters and their sizes.
    set clusterList [getClusters $neighListList]
    set clusterNum [llength $clusterList]
    if {$clusterNum == 0} {
    puts "ERROR: Found zero clusters. Problem with the selection text?"
    exit
    }
    set sizeList {}
    foreach cluster $clusterList {
    lappend sizeList [llength $cluster]
    }
    set sizeList [lsort -integer $sizeList]

    
    set radiusList {}
    set comList {}
    set ligList {}
    
    # code for radius of clusters, what ligands are in each of them
    if 0 {
    foreach cluster $clusterList {

        set sel [atomselect top "index [join $cluster]"]
        set radius [measure rgyr $sel]
        if {$radius > 100} {

            pbc wrap -centersel "index [lindex $cluster 0]" -center com -compound res -now

            set new_radius [measure rgyr $sel]
            set smallest [expr {min($new_radius, $radius)}]

            set com [measure center $sel weight mass]

            lappend radiusList $smallest
            lappend comList $com

            set sphereRadius [expr $smallest * 2]
            set x [lindex $com 0]
            set y [lindex $com 1]
            set z [lindex $com 2]

            set ligandsInCluster [atomselect top "sqrt((x - $x)^2 + (y - $y)^2 + (z - $z)^2) < $sphereRadius"]
            # Get the atom indices
            set indices [$ligandsInCluster get index]

            lappend ligList $indices

            if {$smallest > 300} {
                puts "PROBLEM: $smallest"
                set print "index [join $cluster]"
                puts "cluster: $print"
            }
            pbc unwrap -now
            pbc wrap -compound res -now
        } else {
            set com [measure center $sel weight mass]

            lappend comList $com
            lappend radiusList $radius

            set sphereRadius [expr $radius * 2]
            set x [lindex $com 0]
            set y [lindex $com 1]
            set z [lindex $com 2]

            set ligandsInCluster [atomselect top "sqrt((x - $x)^2 + (y - $y)^2 + (z - $z)^2) < $sphereRadius"]

            # Get the atom indices
            set indices [$ligandsInCluster get index]

            lappend ligList $indices

            if {$radius > 300} {
                    puts "PROBLEM: $radius"
                    set print "index [join $cluster]"
                    puts "cluster: $print"
                }
        }
    }
    }

# Get the largest size
set maxsize 1
set maxclusterList 1
foreach cluster $clusterList {
set size [llength $cluster]
    if { $size > $maxsize} {
    set maxsize $size
    set maxclusterList $cluster
    }
    }

    # Write the sizes.
    puts $outDist $sizeList
    set meanSize [mean $sizeList]

    puts $out $meanSize
    puts $outMax [lindex $sizeList end]
    puts $outN $clusterNum
    puts $outList $clusterList
    puts $outMaxList $maxclusterList
    
    #puts $outComList $comList
    #puts $outRadList $radiusList
    #puts $outLigList $ligList

    # if {$f % $displayPeriod == 0} {
    # puts "FRAME $f: $t clusterNum $clusterNum meanSize $meanSize minSize [lindex $sizeList 0] maxSize [lindex $sizeList end]"
    # }
    
}
# set nFrames0 [expr {$nFrames+$nFrames0}]

close $out
close $outN
close $outMax
close $outDist
close $outMaxList
close $outList
#close $outComList
#close $outRadList
#close $outLigList

# mol delete top
# exit
}
