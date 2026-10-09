#!/usr/bin/env python3
"""
[DEPRECATED] Find HERACLES snapshot closest to target time in days.

DEPRECATED (2026-07-09): HERACLES does not fit the science case for this project
(interacting Type IIn supernovae). The project now uses the STELLA -> CMFGEN
pipeline exclusively. This file is retained for reference only.

For the active snapshot selection path, STELLA .swd files are used directly
via convert_stella_swd_to_cmfgen.py with --target-days.

---

Find HERACLES snapshot closest to target time in days.

Usage:
    python find_closest_snapshot.py <heracles_dir> <target_days>
"""

import sys
import os
import argparse
from pathlib import Path
import re


def find_heracles_snapshots(heracles_dir):
    """
    Find all HERACLES snapshot directories.
    
    Parameters:
    -----------
    heracles_dir : str or Path
        Path to HERACLES output directory
    
    Returns:
    --------
    snapshots : list of dict
        List of snapshot info: {'id': int, 'path': Path, 'time_days': float}
    """
    heracles_dir = Path(heracles_dir)
    snapshots = []
    
    # Look for numeric directories (snapshot IDs)
    for item in heracles_dir.iterdir():
        if item.is_dir() and item.name.isdigit():
            snap_id = int(item.name)
            snap_path = item
            
            # Try to read simulation time from HERACLES output
            # Method 1: Check for .bin file and read time using rd_heracles
            # Method 2: Parse from directory name if it contains time info
            # Method 3: Use directory ID as proxy (not accurate but fallback)
            
            # For now, use directory ID as a proxy
            # In a full implementation, we would read the actual time from the snapshot
            # This requires rd_heracles which needs HERACLES_ROOT setup
            
            # Try to estimate time from directory structure
            # HERACLES typically outputs at regular intervals
            # We'll need to read the actual time from the snapshot data
            
            # For simplicity, we'll return the ID and let the caller handle time extraction
            # In practice, you'd use rd_heracles to get the actual simulation time
            snapshots.append({
                'id': snap_id,
                'path': snap_path,
                'time_days': None  # Will be filled by reading actual data
            })
    
    # Sort by snapshot ID
    snapshots.sort(key=lambda x: x['id'])
    
    return snapshots


def read_snapshot_time(snap_path, heracles_root=None):
    """
    Read actual simulation time from HERACLES snapshot.
    
    Parameters:
    -----------
    snap_path : Path
        Path to snapshot directory
    heracles_root : str or Path, optional
        Path to HERACLES root directory
    
    Returns:
    --------
    time_days : float or None
        Simulation time in days, or None if cannot be read
    """
    try:
        # Try to use rd_heracles to read actual time
        if heracles_root:
            import sys
            visu1 = os.path.join(heracles_root, "heracles-2.0", "visu", "python", "data")
            visu2 = os.path.join(heracles_root, "heracles-2.0", "visu", "python")
            if str(visu1) not in sys.path:
                sys.path.insert(0, str(visu1))
            if str(visu2) not in sys.path:
                sys.path.insert(0, str(visu2))
            
            import rd_heracles as rdh
            snap_id = int(snap_path.name)
            print(f"[DEBUG] Trying to read snapshot {snap_id} from {snap_path.parent}")
            print(f"[DEBUG] snap_path.parent: {snap_path.parent}")
            print(f"[DEBUG] snap_path.parent absolute: {snap_path.parent.absolute()}")
            
            # Try HeraclesData
            try:
                d = rdh.HeraclesData(idump=snap_id, form='bin', verbose=False,
                                     dir=str(snap_path.parent.absolute()))
                # Get time from HERACLES data (correct attribute is 'time' not 't')
                time_s = float(d.time)
                time_days = time_s / 86400.0
                print(f"[DEBUG] Successfully read time with HeraclesData: {time_days:.3f} days")
                return time_days
            except Exception as e:
                print(f"[DEBUG] Failed to read with HeraclesData: {e}")
    except Exception as e:
        print(f"[DEBUG] Exception in read_snapshot_time: {e}")
        # If rd_heracles fails, try alternative methods
        pass
    
    # Fallback: try to parse from output files
    try:
        # Look for .out or .log files with time information
        for file in snap_path.parent.glob(f"heracles_*_{snap_path.name}.out"):
            content = file.read_text()
            # Try to find time information
            match = re.search(r'time\s*=\s*([\d.eE+-]+)\s*s', content)
            if match:
                time_s = float(match.group(1))
                return time_s / 86400.0
    except:
        pass
    
    # If all methods fail, return None
    return None


def find_closest_snapshot(heracles_dir, target_days, heracles_root=None):
    """
    Find HERACLES snapshot closest to target time.
    
    Parameters:
    -----------
    heracles_dir : str or Path
        Path to HERACLES output directory
    target_days : float
        Target time in days
    heracles_root : str or Path, optional
        Path to HERACLES root directory for reading snapshot times
    
    Returns:
    --------
    closest : dict or None
        Closest snapshot info: {'id': int, 'path': Path, 'time_days': float, 'diff_days': float}
    """
    heracles_dir = Path(heracles_dir)
    
    print(f"[DEBUG] heracles_dir: {heracles_dir}")
    print(f"[DEBUG] heracles_dir absolute: {heracles_dir.absolute()}")
    print(f"[DEBUG] heracles_dir exists: {heracles_dir.exists()}")
    print(f"[DEBUG] target_days: {target_days}")
    print(f"[DEBUG] heracles_root: {heracles_root}")
    print(f"[DEBUG] Current working directory: {os.getcwd()}")
    
    if not heracles_dir.exists():
        print(f"ERROR: HERACLES directory not found: {heracles_dir}")
        return None
    
    # Find all snapshots
    snapshots = find_heracles_snapshots(heracles_dir)
    
    if not snapshots:
        print(f"ERROR: No snapshots found in {heracles_dir}")
        return None
    
    print(f"Found {len(snapshots)} snapshots in {heracles_dir}")
    
    # Read actual times from snapshots
    for snap in snapshots:
        print(f"[DEBUG] Processing snapshot {snap['id']}: {snap['path']}")
        snap['time_days'] = read_snapshot_time(snap['path'], heracles_root)
        print(f"[DEBUG] Read time: {snap['time_days']}")
        
        # If time reading failed, skip this snapshot
        if snap['time_days'] is None:
            print(f"  ERROR: Could not read time for snapshot {snap['id']}, skipping")
            continue
    
    # Find closest snapshot
    closest = None
    min_diff = float('inf')
    
    for snap in snapshots:
        if snap['time_days'] is not None:
            diff = abs(snap['time_days'] - target_days)
            print(f"[DEBUG] Snapshot {snap['id']}: {snap['time_days']:.3f} days, diff from target: {diff:.3f} days")
            if diff < min_diff:
                min_diff = diff
                closest = snap
                closest['diff_days'] = diff
    
    if closest:
        print(f"Closest snapshot: {closest['id']} at {closest['time_days']:.2f} days (diff: {closest['diff_days']:.2f} days)")
        print(f"  Path: {closest['path']}")
    else:
        print("ERROR: Could not determine snapshot times")
    
    return closest


def main():
    parser = argparse.ArgumentParser(
        description='Find HERACLES snapshot closest to target time'
    )
    parser.add_argument('heracles_dir', type=str, help='HERACLES output directory')
    parser.add_argument('target_days', type=float, help='Target time in days')
    parser.add_argument('--heracles-root', type=str, default=None,
                       help='HERACLES root directory (for reading snapshot times)')
    parser.add_argument('--output', type=str, default=None,
                       help='Output file to write snapshot path (optional)')
    
    args = parser.parse_args()
    
    # Find closest snapshot
    closest = find_closest_snapshot(args.heracles_dir, args.target_days, args.heracles_root)
    
    if closest:
        # Print snapshot path to stdout for easy capture
        print(str(closest['path']))
        
        # Optionally write to file
        if args.output:
            with open(args.output, 'w') as f:
                f.write(str(closest['path']) + '\n')
                f.write(f"{closest['time_days']:.6f}\n")
            print(f"Snapshot info written to: {args.output}")
        
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
