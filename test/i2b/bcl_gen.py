#!/usr/bin/python2.7

from __future__ import print_function
import random
import os
from datetime import datetime
from xml.etree import ElementTree
from xml.dom import minidom
import struct
import subprocess
import sys
import getopt
import math

#### ONLY EDIT PARAMETERS BETWEEN THESE LINES ####
PARAMS = {
        "lanes": 4,
        "surfaces": 2,
        "swaths": 3,
        "tiles": 12,
        "sections": 3,
        "clusters": 2000,
        "flowcellname": "TESTFLOWCELL",
        "date": datetime.now().strftime("%y%m%d"),
        "reads": [{"num_cycles": 1, "is_indexed": False},
                  {"num_cycles": 1, "is_indexed": False}],
        "dims": {"width": 2048, "height": 7241}
        }
#### ONLY EDIT PARAMETERS BETWEEN THESE LINES ####


# Generates .bcl, .bci, .filter, .locs, RunInfo.xml files to imitate an
# Illumina sequencing machine's output files and directory structure.
# Uses a mix of the bcl2fastq user guide, user guides for each machine,
# and documentation about the .(c)locs file formats from the Picard project
# on GitHub.
#
# author: an8@sanger
#         aacn500@github

machinetypes = [
        "nextseq",
        "hiseqx",
        "hiseq4000",
        "miseq",
        "hiseq2500"
        ]

machinenames = {
        "nextseq": "NSTESTMACHINE",
        "hiseqx": "HXTESTMACHINE",
        "hiseq4000": "HFTESTMACHINE",
        "miseq": "MSTESTMACHINE",
        "hiseq2500": "HSTESTMACHINE",
        }

application_names = {
        "nextseq": "NextSeq Control Software",
        "hiseqx": "HiSeq Control Software",
        "hiseq4000": "HiSeq Control Software",
        "hiseq2500": "HiSeq Control Software",
        "miseq": "MiSeq Control Software"
        }

application_versions = {
        "nextseq": "2.0.0.24",
        "hiseqx": "3.3.39",
        "hiseq4000": "3.3.39",
        "hiseq2500": "2.0.12.0",
        "miseq": "2.5.0.5"
        }

max_params = {
        "nextseq": {
            "lanes": 4,
            "surfaces": 2,
            "swaths": 3,
            "tiles": 12,
            "sections": 3,
            },
        "hiseqx": {
            "lanes": 8,
            "surfaces": 2,
            "swaths": 2,
            "tiles": 24,
            "sections": 1
            },
        "hiseq4000": {
            "lanes": 8,
            "surfaces": 2,
            "swaths": 2,
            "tiles": 28,
            "sections": 1
            },
        "miseq": {
            "lanes": 1,
            "surfaces": 2,
            "swaths": 1,
            "tiles": 19,
            "sections": 1
            },
        "hiseq2500": {
            "lanes": 8,
            "surfaces": 2,
            "swaths": 3,
            "tiles": 16,
            "sections": 1
            }
        }


class Run(object):
    def __init__(self, machinetype):
        assert machinetype in machinetypes
        self.machinetype = machinetype
        self.lanes = [Lane(lane) for lane in xrange(PARAMS["lanes"])]
        self.reads = [Read(read["num_cycles"], read["is_indexed"])
                      for read in PARAMS["reads"]]
        self.id = "{}_{}_{:04d}".format(PARAMS["date"],
                                        machinenames[machinetype],
                                        random.randint(0, 9999))
        self.infopath = ""

    def make_runinfo(self, machinetype):
        root = ElementTree.Element("RunInfo", {
                "xmlns:xsd": "http://www.w3.org/2001/XMLSchema",
                "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
                "Version": "4"
                })

        run = ElementTree.SubElement(root, "Run", {
                "Id": "{}".format(self.id),
                "Number": "2"
                })

        flowcell = ElementTree.SubElement(root, "Flowcell")
        flowcell.text = PARAMS["flowcellname"]

        instrument = ElementTree.SubElement(root, "Instrument")
        instrument.text = machinenames[machinetype]

        date = ElementTree.SubElement(run, "Date")
        date.text = PARAMS["date"]

        reads = ElementTree.SubElement(run, "Reads")
        for i in xrange(len(self.reads)):
            read = ElementTree.SubElement(reads, "Read", {
                    "Number": str(i + 1),
                    "NumCycles": "{:d}".format(self.reads[i].num_cycles),
                    "IsIndexedRead": "Y" if self.reads[i].is_indexed
                                     else "N"
                    })

        flowcellparams = {
            "LaneCount": "{:d}".format(PARAMS["lanes"]),
            "SurfaceCount": "{:d}".format(PARAMS["surfaces"]),
            "SwathCount": "{:d}".format(PARAMS["swaths"]),
            "TileCount": "{:d}".format(PARAMS["tiles"])}

        if machinetype == "nextseq":
            flowcellparams["SectionPerLane"] = "{:d}".format(
                                                        PARAMS["sections"])
            flowcellparams["LanePerSection"] = "2"

        flowcell = ElementTree.SubElement(run, "FlowcellLayout",
                                          flowcellparams)

        if machinetype == "nextseq" or machinetype == "hiseqx" or\
            machinetype == "hiseq4000":
            tile_names = "FiveDigit" if machinetype == "nextseq" \
                                     else "FourDigit"
            tileset = ElementTree.SubElement(flowcell, "TileSet", {
                "TileNamingConvention": tile_names})

            tiles = ElementTree.SubElement(tileset, "Tiles")

            for lane_idx in xrange(PARAMS["lanes"]):
                for section_idx in xrange(PARAMS["sections"]):
                    for swath_idx in xrange(PARAMS["swaths"]):
                        for surface_idx in xrange(PARAMS["surfaces"]):
                            for tile_idx in xrange(PARAMS["tiles"]):
                                tile = ElementTree.SubElement(tiles, "Tile")
                                if machinetype == "nextseq":
                                # Section index is given by number of camera
                                # imaging that section.
                                #
                                #       L1     L2   |   L3     L4
                                #      +---+  +---+ |  +---+  +---+
                                # cam1 |   |  |   | |  |   |  |   | cam4
                                #      |   |  |   | |  |   |  |   |
                                #      |   |  |   | |  |   |  |   |
                                #  ----+---+--+---+-+--+---+--+---+---
                                # cam2 |   |  |   | |  |   |  |   | cam5
                                #      |   |  |   | |  |   |  |   |
                                #      |   |  |   | |  |   |  |   |
                                #  ----+---+--+---+-+--+---+--+---+---
                                # cam3 |   |  |   | |  |   |  |   | cam6
                                #      |   |  |   | |  |   |  |   |
                                #      |   |  |   | |  |   |  |   |
                                #      +---+  +---+ |  +---+  +---+
                                #
                                    section_offset = 1 if lane_idx < 2 else 4
                                    tile.text = "{:d}_{:d}{:d}{:d}{:02d}".\
                                            format(
                                        lane_idx + 1,
                                        surface_idx + 1,
                                        swath_idx + 1,
                                        section_idx + section_offset,
                                        tile_idx + 1)
                                else:
                                    # other machines have no concept of section
                                    tile.text = "{:d}_{:d}{:d}{:02d}".format(
                                            lane_idx + 1,
                                            surface_idx + 1,
                                            swath_idx + 1,
                                            tile_idx + 1)

            image_dims = ElementTree.SubElement(run, "ImageDimensions", {
                "Width": str(PARAMS["dims"]["width"]),
                "Height": str(PARAMS["dims"]["height"])
                })

            image_chans = ElementTree.SubElement(run, "ImageChannels")

            red = ElementTree.SubElement(image_chans, "Name")
            red.text = "Red"

            green = ElementTree.SubElement(image_chans, "Name")
            green.text = "Green"

        if machinetype == "hiseqx" or machinetype == "hiseq2500":
            aligntophix = ElementTree.SubElement(run, "AlignToPhiX")
            if machinetype == "hiseq2500":
                for lane in self.lanes:
                    alignlane = ElementTree.SubElement(aligntophix, "Lane")
                    alignlane.text = "{:d}".format(lane.idx + 1)

        f = open(os.path.join(self.infopath, 'RunInfo.xml'), 'w')
        f.write(pp(root))
        f.close()

    def make_runparams(self, machinetype):
        root = ElementTree.Element("RunParameters", {
            "xmlns:xsd": "http://www.w3.org/2001/XMLSchema",
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance"
            })

        setup = ElementTree.SubElement(root, "Setup")

        name = ElementTree.SubElement(setup, "ApplicationName")
        name.text = application_names[machinetype]

        version = ElementTree.SubElement(setup, "ApplicationVersion")
        version.text = application_versions[machinetype]

        exp_name = ElementTree.SubElement(setup, "ExperimentName")
        exp_name.text = "TestDataExperiment"

        comp_name = ElementTree.SubElement(setup, "ComputerName")
        comp_name.text = "TC1" # test computer 1

        run_date = ElementTree.SubElement(setup, "RunStartDate")
        run_date.text = PARAMS["date"]

        f = open(os.path.join(self.infopath, 'runParameters.xml'), 'w')
        f.write(pp(root))
        f.close()

    def _make_nextseq_bcls(self):
        for lane in self.lanes:
            cycle_idx = 0
            for read in self.reads:
                for cycle in xrange(read.num_cycles):
                    cycle_idx += 1
                    fn = "{:04d}.bcl".format(cycle_idx)
                    f = open(os.path.join(lane.bcpath, fn), mode='wb')
                    f.write(lane.nextseq_bcl())
                    f.close()
                    try:
                        subprocess.call(["bgzip",
                                     os.path.join(lane.bcpath, fn)])
                    except OSError:
                        print("ERROR:",
                            "Couldn't compress bcl file as a blocked gzf.",
                            "Do you have the `bgzip` binary on your PATH?")
                        sys.exit(2)
                    subprocess.call(["mv",
                                     os.path.join(lane.bcpath, fn + '.gz'),
                                     os.path.join(lane.bcpath, fn + '.bgzf')])

    def _make_hiseqx_bcls(self):
        for lane in self.lanes:
            cycle_idx = 0
            for read in self.reads:
                for cycle in xrange(read.num_cycles):
                    cycle_idx += 1
                    for section in lane.sections:
                        for swath in section.swaths:
                            for surface in swath.surfaces:
                                for tile in surface.tiles:
                                    fn = "s_{:d}_{:d}{:d}{:02d}.bcl"\
                                            .format(lane.idx + 1,
                                                    surface.idx + 1,
                                                    swath.idx + 1,
                                                    tile.idx + 1)
                                    f = open(os.path.join(lane.bcpath,
                                             "C{:d}.1".format(cycle_idx),
                                             fn), 'wb')
                                    f.write(tile.hiseqx_bcl())
                                    f.close()
                                    subprocess.call(["gzip",
                                        os.path.join(lane.bcpath,
                                            "C{:d}.1".format(cycle_idx),
                                            fn)])

    def make_bcls(self, machinetype):
        print("Making bcl files...")
        if machinetype == "nextseq":
            self._make_nextseq_bcls()
        elif machinetype == "hiseqx" or machinetype == "hiseq4000":
            self._make_hiseqx_bcls()
        elif machinetype == "miseq":
            self._make_hiseqx_bcls()
        elif machinetype == "hiseq2500":
            self._make_hiseqx_bcls()

    def _make_nextseq_bcis(self):
        # only nextseq generate bci files. One per lane, placed at
        # Data/Intensities/BaseCalls/L00X
        cn = struct.pack("<I", PARAMS["clusters"])
        for lane_idx in xrange(PARAMS["lanes"]):
            lane = self.lanes[lane_idx]
            f = open(os.path.join(lane.bcpath,
                                  's_{:d}.bci'.format(lane_idx + 1)), 'wb')
            s = ''
            for section_idx in xrange(PARAMS["sections"]):
                for swath_idx in xrange(PARAMS["swaths"]):
                    for surface_idx in xrange(PARAMS["surfaces"]):
                        for tile_idx in xrange(PARAMS["tiles"]):
                            section_offset = 1 if lane_idx < 2 else 4
                            s += struct.pack('<I',
                                  int("{:d}{:d}{:d}{:02d}".format(
                                    surface_idx + 1,
                                    swath_idx + 1,
                                    section_idx + section_offset,
                                    tile_idx + 1)))
                            s += cn
            f.write(s)
            f.close()

    def make_bcis(self, machinetype):
        if machinetype == "nextseq":
            print("Making bci files...")
            self._make_nextseq_bcis()
        # other machines don't output bci files
        return

    def _make_nextseq_filters(self):
        # one file per lane, placed at Data/Intensities/BaseCalls/L00X
        for lane_idx in xrange(PARAMS["lanes"]):
            lane = self.lanes[lane_idx]
            s = struct.pack("<I", 0)
            s += struct.pack("<I", 3)
            cluster_passes = ""
            c = 0
            for section in lane.sections:
                for swath in section.swaths:
                    for surface in swath.surfaces:
                        for tile in surface.tiles:
                            for cluster in tile.clusters:
                                c += 1
                                cluster_passes += struct.pack(
                                                    "B", cluster.filtered)
            s += struct.pack("<I", c)
            s += cluster_passes
            f = open(os.path.join(lane.bcpath,
                                  "s_{:d}.filter".format(lane_idx + 1)), "wb")
            f.write(s)
            f.close()

    def _make_hiseqx_filters(self):
        # all machines except nextseq use same filter file format
        # placed at Data/Intensities/BaseCalls/L00X/
        # one file per tile per lane
        for lane_idx in xrange(PARAMS["lanes"]):
            lane = self.lanes[lane_idx]
            for section in lane.sections:
                for swath in section.swaths:
                    for surface in swath.surfaces:
                        for tile in surface.tiles:
                            s = struct.pack("<I", 0)
                            s += struct.pack("<I", 3)
                            cluster_passes = ""
                            c = 0
                            for cluster in tile.clusters:
                                c += 1
                                cluster_passes += struct.pack(
                                                    "B", cluster.filtered)
                            s += struct.pack("<I", c)
                            s += cluster_passes
                            f = open(os.path.join(lane.bcpath,
                                "s_{:d}_{:d}{:d}{:02d}.filter".format(
                                    lane_idx + 1,
                                    surface.idx + 1,
                                    swath.idx + 1,
                                    tile.idx + 1)), 'wb')
                            f.write(s)
                            f.close()

    def make_filters(self, machinetype):
        print("Making filters file...")
        if machinetype == "nextseq":
            self._make_nextseq_filters()
        elif machinetype == "hiseqx" or machinetype == "hiseq4000":
            self._make_hiseqx_filters()
        elif machinetype == "miseq":
            self._make_hiseqx_filters()
        elif machinetype == "hiseq2500":
            self._make_hiseqx_filters()

    def _make_nextseq_locs(self):
        # one locs file per lane, placed in Data/Intensities/L00X/
        for lane_idx in xrange(len(self.lanes)):
            lane = self.lanes[lane_idx]
            s = struct.pack('<I', 1)
            s += struct.pack('<f', 1.0)
            cluster_locs = ""
            c = 0
            for section in lane.sections:
                for swath in section.swaths:
                    for surface in swath.surfaces:
                        for tile in surface.tiles:
                            for cluster in tile.clusters:
                                c += 1
                                cluster_locs += struct.pack(
                                        "<f", random.uniform(0,
                                            PARAMS["dims"]["width"]))
                                cluster_locs += struct.pack(
                                        "<f", random.uniform(0,
                                            PARAMS["dims"]["height"]))
            s += struct.pack("<I", c)
            s += cluster_locs

            f = open(os.path.join(lane.locspath,
                                  "s_{:d}.locs".format(lane_idx + 1)), "wb")
            f.write(s)
            f.close()

    def _make_hiseqx_locs(self):
        # clusters must exist at one of pre-defined "wells", which are in the
        # same locations on each tile, i.e. one locs file for whole run
        # placed in root/Data/Intensities/
        total_clusters = PARAMS["clusters"]
        x_locs = [struct.pack("<f", x) for x in
                        sorted([(random.uniform(0, PARAMS["dims"]["width"]))
                              for xloc in range(total_clusters / 3)])]
        y_locs = [struct.pack("<f", float(y)) for y in
                        sorted([(random.uniform(0, PARAMS["dims"]["height"]))
                              for yloc in range(2 * total_clusters / 3)])]
        cluster_locs = ""
        for xloc in x_locs:
            for yloc in y_locs:
                cluster_locs += xloc
                cluster_locs += yloc

        s = struct.pack('<I', 1)
        s += struct.pack('<f', 1.0)
        s += struct.pack('<I', total_clusters)
        s += cluster_locs

        f = open(os.path.join(self.infopath, 'Data', 'Intensities', 's.locs'),
                 'wb')
        f.write(s)
        f.close()

    def _make_hiseq2500_clocs(self):
        # one clocs file per tile per lane, placed in Data/Intensities/L00X/
        # https://github.com/broadinstitute/picard/blob/master/src/main/java
        #     /picard/illumina/parser/readers/ClocsFileReader.java
        for lane in self.lanes:
            for section in lane.sections:
                for swath in section.swaths:
                    for surface in swath.surfaces:
                        for tile in surface.tiles:
                            # First byte in file gives clocs version (1)
                            s = struct.pack('B', 1)
                            bin_size = 25.0
                            # Should use width given in PARAMS, but
                            # hiseq2500 flowcells have width 2048,
                            # having a different width /might/ break something
                            # (or might not)
                            image_width = 2048 #PARAMS["dims"]["width"]
                            max_x_bins = math.ceil(image_width / bin_size)
                            max_y_bins = 5
                            num_bins = max_x_bins * max_y_bins
                            remaining_clusters = PARAMS["clusters"]

                            # bytes 1-4: unsigned int num_bins
                            s += struct.pack("<I", int(num_bins))
                            for y in xrange(int(num_bins)):
                                cl = min(int(
                                    math.ceil(PARAMS["clusters"] / num_bins)),
                                          remaining_clusters)
                                s += struct.pack("B", cl)
                                remaining_clusters -= cl
                                for cluster in xrange(cl):
                                    s += struct.pack("B", random.randint(
                                                         0, 10 * bin_size - 1))
                                    s += struct.pack("B", random.randint(
                                                         0, 10 * bin_size - 1))
                            f = open(os.path.join(lane.locspath,
                                "s_{:d}_{:d}{:d}{:02d}.clocs".format(
                                                            lane.idx + 1,
                                                            surface.idx + 1,
                                                            swath.idx + 1,
                                                            tile.idx + 1)),
                                                        'wb')
                            f.write(s)
                            f.close()

    def _make_miseq_locs(self):
        # one locs file per tile per lane, placed in Data/Intensities/L00X/
        for lane in self.lanes:
            for section in lane.sections:
                for swath in section.swaths:
                    for surface in swath.surfaces:
                        for tile in surface.tiles:
                            # bytes 0-3: unsigned int locs_version
                            s = struct.pack('<I', 1)
                            # bytes 4-7: float (1.0)
                            s += struct.pack('<f', 1.0)
                            cluster_locs = ""
                            c = 0
                            for cluster in xrange(len(tile.clusters)):
                                c += 1
                                cluster_locs += struct.pack("<f",
                                    random.uniform(0, PARAMS["dims"]["width"]))
                                cluster_locs += struct.pack("<f",
                                    random.uniform(0, PARAMS["dims"]["width"]))
                            # bytes 8-11: unsigned int num_clusters
                            s += struct.pack("<I", c)
                            # bytes 12-end: float x_coord; float y_coord
                            s += cluster_locs

                            f = open(os.path.join(lane.locspath,
                                 "s_{:d}_{:d}{:d}{:02d}.locs".format(
                                                        lane.idx + 1,
                                                        surface.idx + 1,
                                                        swath.idx + 1,
                                                        tile.idx + 1)), 'wb')
                            f.write(s)
                            f.close()

    def make_locs(self, machinetype):
        print("Creating locs file...")
        if machinetype == "nextseq":
            self._make_nextseq_locs()
        elif machinetype == "hiseqx" or machinetype == "hiseq4000":
            self._make_hiseqx_locs()
        elif machinetype == "miseq":
            self._make_miseq_locs()
        elif machinetype == "hiseq2500":
            self._make_hiseq2500_clocs()


class Read(object):
    def __init__(self, num_cycles, is_indexed):
        self.num_cycles = num_cycles
        self.is_indexed = is_indexed


class Lane(object):
    def __init__(self, idx):
        self.sections = [Section(section) for section in xrange(
                                                          PARAMS["sections"])]
        self.bcpath = ""
        self.idx = idx

    def nextseq_bcl(self):
        s = ""
        for section in self.sections:
            for swath in section.swaths:
                for surface in swath.surfaces:
                    for tile in surface.tiles:
                        for cluster in tile.clusters:
                            s += cluster.byte

        l = len(s)
        return struct.pack("<I", l) + s


class Section(object):
    def __init__(self, idx):
        self.swaths = [Swath(swath) for swath in xrange(PARAMS["swaths"])]
        self.idx = idx


class Swath(object):
    def __init__(self, idx):
        self.surfaces = [Surface(surface) for surface in xrange(
                                                          PARAMS["surfaces"])]
        self.idx = idx


class Surface(object):
    def __init__(self, idx):
        self.tiles = [Tile(tile) for tile in xrange(PARAMS["tiles"])]
        self.idx = idx


class Tile(object):
    def __init__(self, idx):
        self.clusters = [Cluster() for cluster in xrange(PARAMS["clusters"])]
        self.idx = idx

    def as_bcl(self):
        for cluster in self.clusters:
            s += cluster.byte
        return s

    def hiseqx_bcl(self):
        s = ""
        for cluster in self.clusters:
            s += cluster.byte
        l = len(s)
        return struct.pack("<I", l) + s


class Cluster(object):
    def __init__(self):
        self.byte = chr(random.randint(0, 255))
        self.filtered = random.randint(0, 1)


def build_directory_structure(run):
    # Illumina output directory name:
    # date (ddmmyy), machinename, four digit id, "_FC"
    os.mkdir(os.path.join(os.getcwd(), run.id + "_FC"))
    os.chdir(os.path.join(os.getcwd(), run.id + "_FC"))
    run.infopath = os.getcwd()
    for l in xrange(len(run.lanes)):
        lane = run.lanes[l]
        lane.idx = l
        os.makedirs(os.path.join(os.getcwd(), 'Data', 'Intensities',
                                 'L{:03d}'.format(l + 1)))
        lane.locspath = os.path.join(os.getcwd(), 'Data', 'Intensities',
                                  'L{:03d}'.format(l + 1))
        os.makedirs(os.path.join(os.getcwd(), 'Data', 'Intensities',
                                  'BaseCalls', 'L{:03d}'.format(l + 1)))
        lane.bcpath = os.path.join(os.getcwd(), 'Data', 'Intensities',
                                  'BaseCalls', 'L{:03d}'.format(l + 1))
        # only nextseq machines do not create cycle directories in lane.bcpath
        if run.machinetype != "nextseq":
            cycle_idx = 0
            for read in run.reads:
                for cycle in xrange(read.num_cycles):
                    cycle_idx += 1
                    os.makedirs(os.path.join(lane.bcpath,
                                             "C{:d}.1".format(cycle_idx)))


def pp(elem):
    """Return a pretty-printed XML string for given element.
    Taken from:
    https://pymotw.com/2/xml/etree/ElementTree/create.html
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


def usage():
    print("Usage:")
    print(" $ ./bcl_gen.py -m <machine type>")
    print("Where machine type is one of:")
    print(" " + ", ".join(machinetypes))
    sys.exit(2)


def main(argv):
    # http://www.diveintopython.net/scripts_and_streams/
    #  command_line_arguments.html"""
    try:
        opts, args = getopt.gnu_getopt(argv, "m:", [])
    except getopt.GetoptError:
        usage()
    for opt, arg in opts:
        if opt in '-m':
            if arg in machinetypes:
                oor_params = []
                for par in max_params[arg].keys():
                    if max_params[arg][par] < PARAMS[par]:
                        print("value of parameter {} is greater ".format(par) +
                              "than machine's normal capabilities.")
                        if (raw_input("Run with max value instead?  ").lower()\
                                not in ['y', 'yes']):
                            print("exiting...")
                            sys.exit(3)
                            return
                        else:
                            PARAMS[par] = max_params[arg][par]
                run = Run(arg)
                build_directory_structure(run)
                run.make_runinfo(run.machinetype)
                run.make_runparams(run.machinetype)
                run.make_bcls(run.machinetype)
                run.make_bcis(run.machinetype)
                run.make_filters(run.machinetype)
                run.make_locs(run.machinetype)
                return
            else:
                usage()
    else:
        usage()
      
if __name__ == "__main__":
    main(sys.argv[1:])
    sys.exit(0)
