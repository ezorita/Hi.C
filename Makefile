SRC_DIR      = src/
C_DIGEST     = re_digest.c
C_HICPARSE   = parse_contacts.c
C_MERGE      = merge_contacts.c
SRC_DIGEST   = $(addprefix $(SRC_DIR), $(C_DIGEST))
SRC_HICPARSE = $(addprefix $(SRC_DIR), $(C_HICPARSE))
SRC_MERGE    = $(addprefix $(SRC_DIR), $(C_MERGE))

FLAGS = -std=c99 -O3
#FLAGS = -std=c99 -g

all: re_digest parse_contacts merge_contacts

re_digest: $(SRC_DIGEST)
	gcc $(FLAGS) $(SRC_DIGEST) -o $@

parse_contacts: $(SRC_HICPARSE)
	gcc $(FLAGS) $(SRC_HICPARSE) -o $@

merge_contacts: $(SRC_MERGE)
	gcc $(FLAGS) $(SRC_MERGE) -o $@

