import Rappture
import sys
import os
import hashlib as h

rap_xml = Rappture.PyXml(sys.argv[1])


def get_interface_data():

    skeleton = rap_xml['input.(skeleton_dropdown).current'].value
    custom = rap_xml['input.(custom_skeleton).current'].value
    kpoints = rap_xml['input.(kpoints).current'].value

    return (skeleton, custom, kpoints)


def get_proj_dir(skeleton, custom, kpoints):

    if (skeleton == "custom"):
        temp_string = custom
    else:
        skeleton_source = os.getenv('OLCAO_DIR') + "/skl/" + skeleton
        temp_string = open(skeleton_source, "r").read()

    temp_string += kpoints

    dir_name = "proj-" + h.md5(temp_string).hexdigest()

    return (dir_name)


def make_proj_dir(dir_name):

    if (not os.path.exists(dir_name)):
        os.mkdir(dir_name)
        existed = False
    else:
        existed = True

    os.chdir(dir_name)

    return (existed)


def make_input(skeleton, custom, kpoints):

    if (skeleton == "custom"):
        with open("olcao.skl", "w") as o:
            o.write(custom)
    else:
        skeleton_source = os.getenv('OLCAO_DIR') + "/skl/" + skeleton
        os.system("cp " + skeleton_source + "  ./olcao.skl")

    os.system("makeinput -kp " + kpoints)


def execute_jobs():

    print ("Starting Optical Properties")
    os.system("uolcao -optc")
    print ("Finished Optical Properties")

    print ("Starting Density of States")
    os.system("uolcao -dos")
    print ("Finished Density of States")

    print ("Starting Bond Order and Q*")
    os.system("uolcao -bond")
    print ("Finished Bond Order and Q*")

    print ("Starting Symmetric Band Structure")
    os.system("uolcao -sybd")
    print ("Finished Symmetric Band Structure")


def gather_tdos_data():

    data = open("gs_dos-fb.t.plot","r").readlines()

    energy = []
    dos = []
    for line in data[1:]:
        values = line.split()
        energy.append(float(values[0]))
        dos.append(float(values[1]))

    tdos = rap_xml['output.curve(tdos)']
    tdos['about.label'] = "Total Density of States"
    tdos['xaxis.label'] = "Energy"
    tdos['yaxis.label'] = "Total Density of States"
    tdos['xaxis.units'] = "eV"
    tdos['yaxis.units'] = "States/[eV cell]"
    tdos['component.xy'] = (energy, dos)


def gather_optc_data():

    data = open("gs_optc-eb.t.plot","r").readlines()

    optc = []
    energy = []
    eps1 = []
    eps2 = []
    elf = []
    for line in data[1:]:
        values = line.split()
        energy.append(float(values[0]))
        eps1.append(float(values[1]))
        eps2.append(float(values[2]))
        elf.append(float(values[3]))

    optc.append(rap_xml['output.curve(eps1)'])
    optc[0]['about.label'] = "Epsilon 1"
    optc[0]['about.group'] = "UV-Vis Dielectric Function"
    optc[0]['xaxis.label'] = "Energy Difference"
    optc[0]['yaxis.label'] = "Eps1, Eps2, ELF"
    optc[0]['xaxis.units'] = "eV"
    optc[0]['yaxis.units'] = "Intensity"
    optc[0]['component.xy'] = (energy, eps1)
    optc.append(rap_xml['output.curve(eps2)'])
    optc[1]['about.label'] = "Epsilon 2"
    optc[1]['about.group'] = "UV-Vis Dielectric Function"
    optc[1]['xaxis.label'] = "Energy Difference"
    optc[1]['yaxis.label'] = "Eps1, Eps2, ELF"
    optc[1]['xaxis.units'] = "eV"
    optc[1]['yaxis.units'] = "Intensity"
    optc[1]['component.xy'] = (energy, eps2)
    optc.append(rap_xml['output.curve(elf)'])
    optc[2]['about.label'] = "Energy Loss Function"
    optc[2]['about.group'] = "UV-Vis Dielectric Function"
    optc[2]['xaxis.label'] = "Energy Difference"
    optc[2]['yaxis.label'] = "Eps1, Eps2, ELF"
    optc[2]['xaxis.units'] = "eV"
    optc[2]['yaxis.units'] = "Intensity"
    optc[2]['component.xy'] = (energy, elf)


def gather_sybd_data():

    data = open("gs_sybd-fb.plot","r").readlines()

    sybd = []  # List of objects. each is a rappture curve.
    sybd_k = []  # x-axis: k-distance on the path
    sybd_energy = []  # List of lists, each is the band data.

    # Compute dimensions
    num_bands = len(data[0].split()) - 3
    num_kpoints = len(data)

    for band in range(num_bands):
        sybd_energy.append([])
        sybd.append([])

    for line in data:
        values = line.split()
        for band in range(num_bands):
            sybd_energy[band].append(float(values[band + 1]))

        # Record the k distance on the path.
        sybd_k.append(float(values[0]))

    for band in range (num_bands):
        sybd[band] = rap_xml['output.curve(sybd[{0}])'.format(band)]
        sybd[band]['about.label'] = "k-dependent Energy Band"
        sybd[band]['about.group'] = "Symmetric Band Structure"
        sybd[band]['xaxis.label'] = "Wave Vector"
        sybd[band]['yaxis.label'] = "Energy"
        sybd[band]['xaxis.units'] = "1/cm"
        sybd[band]['yaxis.units'] = "eV"
        sybd[band]['yaxis.min'] = "-20"
        sybd[band]['yaxis.max'] = "20"
        sybd[band]['component.xy'] = (sybd_k,sybd_energy[band])


def gather_data():
    gather_tdos_data()
    gather_optc_data()
    gather_sybd_data()


# Main subroutine of the script.
def main():
    (skeleton, custom, kpoints) = get_interface_data()

    dir_name = get_proj_dir(skeleton, custom, kpoints)

    if (not make_proj_dir(dir_name)):
        make_input(skeleton, custom, kpoints)

    execute_jobs()

    gather_data()

    # Return to the "rap" directory.
    os.chdir("..")


# Start script execution.
if (__name__ == '__main__'):
    main()


rap_xml.close()
