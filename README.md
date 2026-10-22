# Project-Design-and-Analysis-of-Algorithms
# Lecture Supervision on a Circular Timetable

Selects a subset S of lectures such that every lecture overlaps exactly one
lecture of S (a lecture never counts itself). Lectures repeat every 24 hours
and may wrap past midnight. Algorithm: O(n^3) time, O(n^2) bits of space
(see `report.pdf`).

## Requirements
- Python 3.8 or newer (standard library only, nothing to install)
- To rebuild the report: a LaTeX distribution (pdflatex, latexmk)

## Files
- `lecture_solver.py` - the algorithm, brute-force checker, command-line tool
- `sample.txt`        - example input
- `report.pdf`        - the submitted project report

## Input format
First line: n. Then n lines, one per lecture: `start end` as `HH:MM` (24-hour).
Text after `#` is ignored.

    4
    23:00 01:00     # lecture 1 wraps past midnight
    00:30 02:00
    13:00 15:00
    14:00 16:00

Conventions:
- Intervals are half-open: a lecture ending at 10:00 and one starting at 10:00
  do NOT overlap.
- `start end` equal (e.g. `00:00 00:00`) means a full 24-hour lecture.
- Lectures are numbered 1..n in the order given.

## Running
    python3 lecture_solver.py                 # type the lectures in interactively
    python3 lecture_solver.py mydata.txt      # read your own file
    cat mydata.txt | python3 lecture_solver.py -   # read from standard input
    python3 lecture_solver.py --demo          # 8 built-in test cases
    python3 lecture_solver.py --stress 1000   # compare with brute force

## Output
Either `No valid subset exists.` or the selected lectures (1-indexed), which
selected lecture supervises each lecture, and an independent verification line
(`PASSED`) that re-checks the answer from scratch.

## GitHub Actions
The workflow in `.github/workflows/test_solver.yml` runs automatically on pushes to `main` and can also be started manually from the repository's Actions tab. It runs the three demonstration cases, the sample input, `mydata.txt`, and 100 randomized stress cases.

## Report
Open `report.pdf` for the submitted project report. The LaTeX source file is not included in this repository.
