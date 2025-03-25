source cluster_size.tcl
package require topotools
package require pbctools

# Define the base path
set base_path "/n/netscratch/shakhnovich_lab/Lab/jiojeong/simulations/20250319_metabolon_subtract_kT_from_COQ_PPI"

# List all subdirectories in the base path
set subdirs [glob -directory $base_path -type d *]

set all_dirs [list]

foreach subdir $subdirs {
    set temp [glob -directory $subdir -type d *]
    foreach t $temp {
        lappend all_dirs $t
    }
}


# Loop through each directory
foreach dir $all_dirs {
    # Define the full path to a file, assuming a known structure or naming
    set system_data_file_path "$dir/system.data"
    set dcd_files [glob -nocomplain -directory $dir equi1*.dcd]
    if {[llength $dcd_files] > 0} {
        set dcd_file_path [lindex $dcd_files 0]}
    # Check if the file exists before trying to process it
    if {[file exists $system_data_file_path]} {
        puts "Processing file: $dir"
        topo readlammpsdata $system_data_file_path 
        #animate read dcd $dcd_file_path beg 0 end -1 skip 50 waitfor all
        animate read dcd $dcd_file_path beg 0 end -1 skip 200 waitfor all
        #animate read dcd $dcd_file_path beg 0 end -1 waitfor all
        pbc wrap -compound res -all
        file mkdir $dir/cluster_results_skip_first_frame
        cluster_size top 
        mol delete top}
}

# Notify completion
puts "All directories processed."
