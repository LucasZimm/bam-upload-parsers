# bam-upload-parsers

Some parsers for openBIS to be used in the bam-masterdata parser interface

This repository is intended to be used as entry point for the parsers from [`openbis-upload-helper`](https://github.com/BAMresearch/openbis-upload-helper).


## Table of Contents

- [bam-upload-parsers](#bam-upload-parsers)
  - [Table of Contents](#table-of-contents)
  - [About](#about)
  - [Features](#features)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation](#installation)

---

## About

This repository restructures parser components originally from `bam-upload-helper`, aiming to improve modularity and maintainability for integration with the bam-masterdata parser interface.

## Features

- parser implementations
- Built for integration with openBIS environments
- Includes automated test suite

## Getting Started

### Prerequisites

- Tested for Python 3.10+
- `pip` for dependency installation
- Optional: virtualenv or conda environment for isolation

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/LucasZimm/bam-upload-parsers.git
   cd bam-upload-parsers
2. Create and activate a virtual environment:
    ```bash
    python -m .venv venv
    source venv/bin/activate  # Windows:    venv\Scripts\activate

3. Install dependencies:
    ```bash
    pip install .

