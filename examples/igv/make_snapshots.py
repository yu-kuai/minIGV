#!/usr/bin/env python3
"""Recreate the example IGV snapshots in ../raw from public ENCODE data.

1. Start IGV desktop (tested with IGV 2.16) and leave it open.
   Its batch port must be enabled: View > Preferences > Advanced >
   "Enable port" (default 60151, on by default).
2. Run:  python make_snapshots.py

IGV streams the ENCODE K562 tracks listed in encode_k562_session.xml from
encodeproject.org (loading takes a minute or two), then saves one SVG per
locus into ../raw.
"""
import os
import socket

HERE = os.path.dirname(os.path.abspath(__file__))
SESSION = os.path.join(HERE, "encode_k562_session.xml")
OUT = os.path.normpath(os.path.join(HERE, "..", "raw"))

LOCI = {
    "HOXA_K562.svg": "chr7:27,090,000-27,250,000",    # Polycomb-repressed HOXA cluster
    "MYC_K562.svg": "chr8:127,730,000-127,745,000",   # highly expressed MYC
    "GAPDH_K562.svg": "chr12:6,530,000-6,545,000",    # housekeeping GAPDH
}


def main(port=60151):
    os.makedirs(OUT, exist_ok=True)
    with socket.create_connection(("127.0.0.1", port), timeout=900) as s:
        f = s.makefile("rw")

        def cmd(c):
            f.write(c + "\n")
            f.flush()
            reply = f.readline().strip()
            print("%-60s %s" % (c[:60], reply))

        cmd("new")
        cmd("load " + SESSION)
        cmd("snapshotDirectory " + OUT)
        cmd("maxPanelHeight 1000")
        for name, locus in LOCI.items():
            cmd("goto " + locus)
            cmd("snapshot " + name)


if __name__ == "__main__":
    main()
