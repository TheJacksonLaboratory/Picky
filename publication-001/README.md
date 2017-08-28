# Publication

Gong et al. "Nanopore Sequencing Reveals High-Resolution Structural Variation in the Cancer Genome" (Submitted)

Each MinION sequenced read is assigned a unique 128-bits number (aka read_uuid) by the MinKNOW software. The read_uuid is reported as 32 hexidecimal digits (e.g. 27bc311c-87ef-496b-9153-7f198ed69e63). 
To make tracking reads alignment easier in IGV browser, we assigned a human-friendly read id composed of 7 elements as follow:

| Element | Description | Examples |
|:---:| --- | --- |
| 1 | experiment id | WTD09 |
| 2 | unique running serial number from 1 assigned to each read | 1 |
| 3 | 1D or 2D status of read | 2D |
| 4 | Pass or Fail status of read | P |
| 5 | Read code | 1 = Read template sequence is present<br>2 = Read complement sequence is present<br>3 = Both read template and complement sequence are present |
| 6 | Read sequence code | T=Template<br>C=Complement<br>M=Merged consensus 2D sequence |
| 7 | Read sequence length prefixed with 'L' | L23284 |

To facilitate tracking these IDs back to the read_uuid, each run/experiment has an associated hash table (\<experiment\>-hash.xls).
The first line of the hash table contains the columns header prefixed by '#'.
Each line contains 7 columns as follow:

| Column | Description | Examples |
|:---:| --- | --- |
| 1 | read_uuid | UUID assigned by MinKNOW software |
| 2 | fast5 | .fast5 file where read sequences is extracted from |
| 3 | PassFail | P=Pass, F=Fail |
| 4 | 1D_2D | 1D read(s) or 2D-reads |
| 5 | template_hash_id | human-friendly read id for read's template sequence |
| 6 | complement_hash_id | human-friendly read id for read's template sequence |
| 7 | 2d_hash_id | human-friendly read id for read's template sequence |

## Experiment Hash Table

| Experiment | Run Type | Hash Table | File Size | MD5 |
|:---:|:---:|:---:|:---:| --- |
| WTD01 | 2D | [WTD01-hash.xls.gz](WTD01-hash.xls.gz) | 7.3M | 10135aeaf5b6cefdd8e6a8b16ef722dc |
| WTD02 | 2D | [WTD02-hash.xls.gz](WTD02-hash.xls.gz) | 5.5M | 24d0f4a1a3bec579c2ddc12b5a6de40f |
| WTD03 | 2D | [WTD03-hash.xls.gz](WTD03-hash.xls.gz) | 2.9M | 6ad27eedc54c061432860d511643bae6 |
| WTD04 | 2D | [WTD04-hash.xls.gz](WTD04-hash.xls.gz) | 8.2M | d5de63b8e452ab705d5bd3e365b0c298 |
| WTD05 | 2D | [WTD05-hash.xls.gz](WTD05-hash.xls.gz) | 13M | 44294abd2715ae58cfb771d22ff01b49 |
| WTD06 | 2D | [WTD06-hash.xls.gz](WTD06-hash.xls.gz) | 9.1M | ea150dff613c7be1d3c0123ec8710628 |
| WTD07 | 2D | [WTD07-hash.xls.gz](WTD07-hash.xls.gz) | 3.6M | 946aa2a2dce7cf30d84dbb6d594f2cd6 |
| WTD08 | 2D | [WTD08-hash.xls.gz](WTD08-hash.xls.gz) | 3.2M | 9953493fca33d13d6d77677be0a3397c |
| WTD09 | 2D | [WTD09-hash.xls.gz](WTD09-hash.xls.gz) | 2.2M | 78b9143b0151040000ae6744969faac1 |
| WTD10 | 2D | [WTD10-hash.xls.gz](WTD10-hash.xls.gz) | 2.6M | aced30aa357832170eb420e9ff97ce47 |
| WTD11 | 2D | [WTD11-hash.xls.gz](WTD11-hash.xls.gz) | 3.9M | 8b0400096dcddc9c2d8a9f8ce0404580 |
| WTD12 | 2D | [WTD12-hash.xls.gz](WTD12-hash.xls.gz) | 4.9M | b30315745e5fec1331b3d6920e79f42f |
| WTD13 | 2D | [WTD13-hash.xls.gz](WTD13-hash.xls.gz) | 7.1M | aee3cf4861d09899b6b59336facf7e0f |
| WTD14 | 1D | [WTD14-hash.xls.gz](WTD14-hash.xls.gz) | 7.1M | 4b9e79adb49e4fdc812716d367884155 |
| WTD15 | 1D | [WTD15-hash.xls.gz](WTD15-hash.xls.gz) | 5.7M | 12b67e2a9899055a32754fdaeae42416 |

All MD5 checksums are also available in a single file, [hash-checksums.md5](hash-checksums.md5).



