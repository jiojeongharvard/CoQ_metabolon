import re
import random
import math
import sys


def compute_distance(coord1, coord2):
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    return math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def generate_random_coordinate(max_coord):
    x1 = random.uniform(0, max_coord)
    y1 = random.uniform(0, max_coord)
    z1 = random.uniform(0, max_coord)
    return (x1, y1, z1)

def generate_coordinates_with_distance(rad, max_coord):
   while True:
        # Generate the first random coordinate as a float
        x1 = random.uniform(0, max_coord)
        y1 = random.uniform(0, max_coord)
        z1 = random.uniform(0, max_coord)

        # Calculate a random direction for the distance vector
        theta = random.uniform(0, 2 * math.pi)
        phi = random.uniform(0, math.pi)
        
        # Calculate the distance vector components
        dx = rad * math.sin(phi) * math.cos(theta)
        dy = rad * math.sin(phi) * math.sin(theta)
        dz = rad * math.cos(phi)

        # Calculate the second coordinate
        x2 = x1 + dx
        y2 = y1 + dy
        z2 = z1 + dz

        # Check if the second coordinate is within bounds. Second condition is to ensure a pair with a distance like 10.0001 is not returned from rounding
        if (0 <= x2 <= max_coord and 0 <= y2 <= max_coord and 0 <= z2 <= max_coord) and (compute_distance((x1, y1, z1), (x2, y2, z2)) <= rad):
            return (x1, y1, z1), (x2, y2, z2)

def is_valid_point(new_point, points, min_distance):
    return all(compute_distance(new_point, point) > min_distance for point in points)

       
if __name__ == "__main__":
    #random.seed(42)
           
    enzArr = []
    num_crowd = 0
    num_lig = 0
    rad_enz = 0
    box_len = 0
    enzMass = 0
    activeSiteMass = 0
    ligMass = 0
    crowdMass = 0
    num_coq9 = 0
    coq9Mass = 0
    ligType = 0
    cutoff_enz_enz = 21
    cutoff_enz_lig = 12
    cutoff_lig_lig = 4
    
    user_input_enzymes = input("Enter distribution of enzyme molecules (ex. 100 100): ")
    match = re.match(r'(\d+(?:\s+\d+)*)', user_input_enzymes)
    if match:
        enzStr = match.group(1)
        enzArr = list(map(int, enzStr.split()))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)

    user_input_enz_mass = input("Enter mass of enzymes (ex. 30000): ")
    match = re.match(r'(^\d+$)', user_input_enz_mass)
    if match:
        enzMass = int(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)
        
    user_input_AS_mass = input("Enter mass of enzyme active site (ex. 1500): ")
    match = re.match(r'(^\d+$)', user_input_AS_mass)
    if match:
        activeSiteMass = int(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)
        
    user_input_lig = input("Enter number of starting ligand molecules (ex. 1000): ")
    match = re.match(r'(^\d+$)', user_input_lig)
    if match:
        num_lig = int(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)

    user_input_lig_mass = input("Enter mass of ligand (ex. 800): ")
    match = re.match(r'(^\d+$)', user_input_lig_mass)
    if match:
        ligMass = int(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)
    
    user_input_lig_type = input("Enter number of unique ligands in entire pathway(ex. 8): ")
    match = re.match(r'(^\d+$)', user_input_lig_type)
    if match:
        ligType = int(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)
        
    user_input_crowder = input("Enter number of crowder molecules (ex. 700): ")
    match = re.match(r'(^\d+$)', user_input_crowder)
    if match:
        num_crowd = int(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)

    user_input_crowd_mass = input("Enter mass of crowder (ex. 30000): ")
    match = re.match(r'(^\d+$)', user_input_crowd_mass)
    if match:
        crowdMass = int(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)
    
    user_input_coq9 = input("Enter number of COQ9 molecules (ex. 50): ")
    match = re.match(r'(^\d+$)', user_input_coq9)
    if match:
        num_coq9 = int(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)

    user_input_coq9_mass = input("Enter mass of COQ9 (ex. 30000): ")
    match = re.match(r'(^\d+$)', user_input_coq9_mass)
    if match:
        coq9Mass = int(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)
        
        
    user_input_radius = input("Enter enzyme radius (ex. 10): ")
    match = re.match(r'(^\d+$)', user_input_radius)
    if match:
        rad_enz = float(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)
        
    user_input_box_length = input("Enter box size (ex. 800): ")
    match = re.match(r'(^\d+$)', user_input_box_length)
    if match:
        box_len = float(match.group(1))
    else:
        print("The input string does not match the expected pattern.")
        sys.exit(1)


    enzType = len(enzArr)
    totalEnzymeMolecules = sum(enzArr)
    totalLigandMolecules = num_lig
    totalMolecules = totalEnzymeMolecules + totalLigandMolecules + num_crowd + num_coq9
    totalAtoms = totalEnzymeMolecules * 2 + totalLigandMolecules + num_crowd + num_coq9
        
    print(f"Total number of Molecules: {totalMolecules}, Total number of Atoms: {totalAtoms}\nDistribution of enzymes: {enzArr}, Number of enzyme types: {enzType}, Mass of enzyme: {enzMass}, Mass of active site: {activeSiteMass}\nNumber of starting ligand: {num_lig}, Mass of ligands: {ligMass}, Number of ligand types in pathway: {ligType}\nNumber of crowder: {num_crowd}, Mass of crowder: {crowdMass}\nNumber of COQ9: {num_coq9}, Mass of COQ9: {coq9Mass} Radius of enzyme: {rad_enz},  Box length: {box_len}")

    totalInt = 2* enzType + ligType + 2
    enzIndex = []
    activeSiteIndex = []
    ligIndex = []
    crowderIndex = []
    coq9Index = []
    for i in range(1, totalInt+1):
        if i <= enzType:
            enzIndex.append(i)
        elif i <= enzType * 2:
            activeSiteIndex.append(i)
        elif i <= ligType+ enzType * 2:
            ligIndex.append(i)
        elif i <= ligType+ enzType * 2 + 1:
            crowderIndex.append(i)
        else:
            coq9Index.append(i)

    print(enzIndex)
    print(activeSiteIndex)
    print(ligIndex)
    print(crowderIndex)
    print(coq9Index)
    
    enz_crowd_coordinates = []
    active_site_coordinates = []
    lig_coordinates = []


    with open('new_system.data', 'w') as file:
        original_stdout = sys.stdout
        sys.stdout = file
        print("LAMMPS Description")
        print()
        print(f"{totalAtoms}  atoms")
        print("0  bonds")
        print("0  angles")
        print("0  dihedrals")
        print("0  impropers")
        print()
        print(f"{totalInt}  atom types")
        print("1 bond types")
        print("20 extra bond per atom")
        print(f"0 {box_len} xlo xhi")
        print(f"0 {box_len} ylo yhi")
        print(f"0 {box_len} zlo zhi")
        print()
        print("Masses")
        print()
        for i in range(len(enzIndex)):
            print(f"{enzIndex[i]} {enzMass}")
        
        for i in range(len(activeSiteIndex)):
            print(f"{activeSiteIndex[i]} {activeSiteMass}")
            
        for i in range(len(ligIndex)):
            print(f"{ligIndex[i]} {ligMass}")
        
        for i in range(len(crowderIndex)):
            print(f"{crowderIndex[i]} {crowdMass}")
        
        for i in range(len(coq9Index)):
            print(f"{coq9Index[i]} {coq9Mass}")
        print()
        print("Atoms")
        print()
        
        atomCounter = 1
        moleculeCounter = 1
        while (moleculeCounter <= totalMolecules):
            if moleculeCounter <= totalEnzymeMolecules:
                for i in range(enzType):
                    curr_enzyme = enzIndex[i]
                    curr_activeSite = activeSiteIndex[i]
                    
                    curr_enzyme_molecule_num = enzArr[i]
                    for j in range(curr_enzyme_molecule_num):
                        while True:
                            coord1, coord2 = generate_coordinates_with_distance(rad_enz, box_len)
                            if is_valid_point(coord1, enz_crowd_coordinates, cutoff_enz_enz):
                                print(f"{atomCounter} {moleculeCounter} {curr_enzyme} 0 {coord1[0]} {coord1[1]} {coord1[2]}")
                                atomCounter += 1
                                print(f"{atomCounter} {moleculeCounter} {curr_activeSite} 0 {coord2[0]} {coord2[1]} {coord2[2]}")
                                atomCounter += 1
                                moleculeCounter += 1
                                enz_crowd_coordinates.append(coord1)
                                active_site_coordinates.append(coord2)
                                break
            
            elif moleculeCounter <= totalEnzymeMolecules + num_coq9:
                for j in range(num_coq9):
                    while True:
                        coord = generate_random_coordinate(box_len)
                        coq9 = coq9Index[0]
                        if is_valid_point(coord, enz_crowd_coordinates, cutoff_enz_enz):
                            print(f"{atomCounter} {moleculeCounter} {coq9} 0 {coord[0]} {coord[1]} {coord[2]}")
                            atomCounter += 1
                            moleculeCounter += 1
                            enz_crowd_coordinates.append(coord)
                            break
                    
            elif moleculeCounter <= totalEnzymeMolecules + num_crowd + num_coq9:
                for j in range(num_crowd):
                    while True:
                        coord = generate_random_coordinate(box_len)
                        crowder = crowderIndex[0]
                        if is_valid_point(coord, enz_crowd_coordinates, cutoff_enz_enz):
                            print(f"{atomCounter} {moleculeCounter} {crowder} 0 {coord[0]} {coord[1]} {coord[2]}")
                            atomCounter += 1
                            moleculeCounter += 1
                            enz_crowd_coordinates.append(coord)
                            break
                    
            else: 
                for i in range(num_lig):
                    while True:
                        curr_lig = ligIndex[0]
                        coord = generate_random_coordinate(box_len)
                        if is_valid_point(coord, enz_crowd_coordinates, cutoff_enz_lig) and is_valid_point(coord, lig_coordinates, cutoff_lig_lig):
                            print(f"{atomCounter} {moleculeCounter} {curr_lig} 0 {coord[0]} {coord[1]} {coord[2]}")
                            atomCounter += 1
                            moleculeCounter += 1
                            lig_coordinates.append(coord)
                            break
                        
        sys.stdout = original_stdout
