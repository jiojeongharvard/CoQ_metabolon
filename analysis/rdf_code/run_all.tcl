package require topotools
package require pbctools

# Define the base path
set base_path "/n/netscratch/shakhnovich_lab/Lab/jiojeong/simulations/20241122_metabolon_CoQ9_HSHS1KT_self_interaction_1KT_vary_kcat_new_dump"

# List all subdirectories in the base path
set subdirs [glob -directory $base_path -type d *]

set all_dirs [list]
set all_dirs2 [list]
puts $subdirs

foreach subdir $subdirs {
    set temp [glob -directory $subdir -type d *]
    foreach t $temp {
        lappend all_dirs $t
    }
}

foreach subdir $all_dirs {
    set temp [glob -directory $subdir -type d *]
    foreach t $temp {
        lappend all_dirs2 $t
    }
}

# Loop through each directory
foreach dir $all_dirs2 {
    # Define paths for system data and DCD files
    set system_data_file_path "$dir/system.data"
    set dcd_files [glob -nocomplain -directory $dir equi1*.dcd]

    # Check if system.data and at least one DCD file exist
    if {[file exists $system_data_file_path] && [llength $dcd_files] > 0} {
        puts "Processing directory: $dir"

        # Load system data
        topo readlammpsdata $system_data_file_path 

        # Process the first available DCD file
        set dcd_file_path [lindex $dcd_files 0]
        puts "Loading trajectory: $dcd_file_path"
        animate read dcd $dcd_file_path beg 0 end -1 skip 200 waitfor all

        # Apply periodic boundary conditions
        pbc wrap -compound res -all

        # Define atom selections
        set sel1_6 [atomselect top "type 1"]  ;# COQ6
        set sel2_3 [atomselect top "type 2"]  ;# COQ3
        set sel3_4 [atomselect top "type 3"]  ;# COQ4
        set sel4_5 [atomselect top "type 4"]  ;# COQ5
        set sel5_7 [atomselect top "type 5"]  ;# COQ7

        # Get simulation box dimensions (assumes an orthorhombic box)
        set dims [molinfo top get {a b c}]
        set V [expr {[lindex $dims 0] * [lindex $dims 1] * [lindex $dims 2]}]

        # Compute densities (atoms per unit volume) for each selection.
        set density_6 [expr {[$sel1_6 num] / double($V)}]
        set density_3 [expr {[$sel2_3 num] / double($V)}]
        set density_4 [expr {[$sel3_4 num] / double($V)}]
        set density_5 [expr {[$sel4_5 num] / double($V)}]
        set density_7 [expr {[$sel5_7 num] / double($V)}]

        # Define parameters for RDF calculation
        set delta_r 0.5   ;# Bin size for RDF
        set r_max 100.0    ;# Maximum radius for RDF

        # Compute time-averaged raw RDF using built-in frame range options:
        # The following commands average over frames from 0 to -1 (last frame) with a step of 100.
        set res_6_3 [measure gofr $sel1_6 $sel2_3 delta $delta_r rmax $r_max first 0 last -1 step 1 usepbc True]
        set res_3_4 [measure gofr $sel2_3 $sel3_4 delta $delta_r rmax $r_max first 0 last -1 step 1 usepbc True]
        set res_4_6 [measure gofr $sel3_4 $sel1_6 delta $delta_r rmax $r_max first 0 last -1 step 1 usepbc True]
        set res_6_5 [measure gofr $sel1_6 $sel4_5 delta $delta_r rmax $r_max first 0 last -1 step 1 usepbc True]
        set res_5_7 [measure gofr $sel4_5 $sel5_7 delta $delta_r rmax $r_max first 0 last -1 step 1 usepbc True]
        set res_7_3 [measure gofr $sel5_7 $sel2_3 delta $delta_r rmax $r_max first 0 last -1 step 1 usepbc True]

        # Each "measure gofr" returns a two-element list: {r_values raw_counts}.
        # We assume r_values is the same for all pairs.
        set r_vals [lindex $res_6_3 0]
        set rdf_6_3 [lindex $res_6_3 1]
        set rdf_3_4 [lindex $res_3_4 1]
        set rdf_4_6 [lindex $res_4_6 1]
        set rdf_6_5 [lindex $res_6_5 1]
        set rdf_5_7 [lindex $res_5_7 1]
        set rdf_7_3 [lindex $res_7_3 1]
        set list_timeframes [lindex $res_6_3 4]
        set num_timeframes [lindex $list_timeframes 0]

        # Determine number of bins
        set num_bins [llength $r_vals]

        # Open the output file and write header
        set outfile [open "$dir/rdf_avg.dat" w]
        puts $outfile "# r   RDF_COQ6-3   RDF_COQ3-4   RDF_COQ4-6   RDF_COQ6-5   RDF_COQ5-7   RDF_COQ7-3"

        # Process each bin: normalize by the ideal gas expectation (shell volume * partner density)
        for {set i 0} {$i < $num_bins} {incr i} {
            set r [lindex $r_vals $i]
            # Spherical shell volume: 4πr²Δr
            set shell_vol [expr {4.0 * 3.141592653589793 * $r * $r * $delta_r}]
            set norm_6_3 [expr { [lindex $rdf_6_3 $i] }]
            set norm_3_4 [expr { [lindex $rdf_3_4 $i] }]
            set norm_4_6 [expr { [lindex $rdf_4_6 $i] }]
            set norm_6_5 [expr { [lindex $rdf_6_5 $i] }]
            set norm_5_7 [expr { [lindex $rdf_5_7 $i] }]
            set norm_7_3 [expr { [lindex $rdf_7_3 $i] }]
            puts $outfile "[format "%.3f" $r] [format "%.5f %.5f %.5f %.5f %.5f %.5f" $norm_6_3 $norm_3_4 $norm_4_6 $norm_6_5 $norm_5_7 $norm_7_3]"
        }

        close $outfile
        puts "Normalized RDF saved to: $dir/rdf_avg.dat"

        # Cleanup: delete the loaded molecule
        mol delete top
    } else {
        puts "Skipping directory: $dir (Missing system.data or DCD file)"
    }
}

puts "All directories processed."
