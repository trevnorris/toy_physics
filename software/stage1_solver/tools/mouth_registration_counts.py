"""Sufficient finite-occupancy test for a snapshot of the renewing contact law."""
import json
import math
from fractions import Fraction
from pathlib import Path

from scipy.stats import poisson


def check():
    archive = Path("research/pde_ledger_v3/_measurements/mouth_collar_charge_transmission.json")
    lines = json.loads(archive.read_text())["runs"]["mouth_collar_attachment_window01"]["stdout_lines"]
    def field(line, key):
        return Fraction(line.split(key+"=")[1].split()[0])
    weak = min(field(line,"transmission_for_g_below") for line in lines)
    strong = max(field(line,"screening_for_g_above") for line in lines)
    level = max(int(field(line,"level")) for line in lines)
    rim_links = 12 * 4**level
    feed_links = 6 * (rim_links + 4*4**level + 2)
    print(f"mesh_level={level} rim_links={rim_links} feed_links={feed_links}")
    for occupancy in (16, 1024, 16384):
        # Each realization has the same mean elastic strengths. Occupancy and
        # single-contact stiffness are separate material inputs. nu=j*tau_feed.
        one_contact = Fraction(2,occupancy)
        rim_minimum = math.ceil(64/one_contact)
        for label in ("delayed","immediate"):
            if label=="delayed":
                rim_mean = occupancy * 128/3
                feed_mean = occupancy/257
                # Need EVERY feed contact group strictly below the weak ceiling.
                first_bad = math.ceil(weak/one_contact)
                feed_bad = poisson.sf(first_bad-1,feed_mean)
            else:
                rim_mean = 128*occupancy
                feed_mean = occupancy
                # EVERY feed link above the strong floor suffices.
                first_good = math.floor(strong/one_contact)+1
                feed_bad = poisson.cdf(first_good-1,feed_mean)
            rim_bad = poisson.cdf(rim_minimum-1,rim_mean)
            # Union bound uses only marginal Poisson counts on the links.
            fail_upper = min(1.,rim_links*rim_bad + feed_links*feed_bad)
            print(f"nu={occupancy} law={label} single_contact_over_K={float(one_contact):.12g} "
                  f"feed_mean_registered={feed_mean:.12g} rim_mean_registered={rim_mean:.12g} "
                  f"feed_bad_per_link={feed_bad:.10g} "
                  f"snapshot_certificate_failure_upper={fail_upper:.10g}")
            assert 0 <= fail_upper <= 1


if __name__ == "__main__":
    check()
