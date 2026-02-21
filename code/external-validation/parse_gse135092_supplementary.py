#!/usr/bin/env python3
"""
Parse Excel supplementary files to find AMD subtype metadata
"""
import xml.etree.ElementTree as ET
import zipfile
import sys
import os

def parse_excel_structure(xlsx_path):
    """Extract structure and sample data from Excel file"""
    print(f"\n{'='*70}")
    print(f"FILE: {os.path.basename(xlsx_path)}")
    print('='*70)
    
    with zipfile.ZipFile(xlsx_path, 'r') as z:
        # Get shared strings
        try:
            shared_xml = z.read('xl/sharedStrings.xml')
            root = ET.fromstring(shared_xml)
            strings = []
            for si in root.findall('.//{http://schemas.openxmlformats.org/spreadsheetml/2006/main}t'):
                if si.text:
                    strings.append(si.text)
            
            print(f"\nTotal unique strings: {len(strings)}")
            
            # Look for AMD-related keywords
            amd_keywords = ['AMD', 'amd', 'control', 'Control', 'donor', 'Donor', 
                           'sample', 'Sample', 'subtype', 'stage', 'geographic', 
                           'neovascular', 'intermediate', 'wet', 'dry', 'atrophy']
            
            found_keywords = []
            for kw in amd_keywords:
                matches = [s for s in strings if kw in s]
                if matches:
                    found_keywords.append((kw, len(matches)))
                    if len(matches) <= 20:
                        print(f"\n  Keyword '{kw}' found in: {matches[:10]}")
            
            if found_keywords:
                print(f"\n  Summary of keywords found:")
                for kw, count in found_keywords:
                    print(f"    {kw}: {count} occurrences")
            else:
                print("\n  No AMD-related keywords found")
                
        except Exception as e:
            print(f"  Error reading shared strings: {e}")

if __name__ == "__main__":
    supp_dir = "/home/ysuhail/work/Tannin-AMD/data/external/geo/GSE135092/supplementary"
    
    for i in range(2, 6):
        xlsx_file = os.path.join(supp_dir, f"mmc{i}.xlsx")
        if os.path.exists(xlsx_file):
            parse_excel_structure(xlsx_file)
