
parser <- argparse::ArgumentParser()

parser $ add_argument("-v", "--verbose", action = "store_true", default = TRUE,
                      help = "print detailed output [default]")
parser $ add_argument("-q", "--quiet", action = "store_false",
                      dest = "verbose", help = "print only necessary outputs.")

parser $ add_argument("-t", "--threads", type = "integer", default = 40,
                      help = "number of parallel workers.", dest = "threads")

# parser $ add_argument(type = "character", dest = "inputs", nargs = "+",
#                       help = "input accessions")

args <- parser $ parse_args()
