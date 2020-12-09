import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

  

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """

    transfer_table = {
                    0: {    0: 1 - PROBS["mutation"],
                            1: 0 + PROBS["mutation"]
                        },
                    1: {    0: 0.5,
                            1: 0.5
                        },
                    2: {    0: 0 + PROBS["mutation"],
                            1: 1 - PROBS["mutation"]
                        }
                }
        
        
    #calculate 1 gene list
    work_dict_first = {}
    work_dict_second = {}
    for person in people:
        person_gene = 0
        person_trait = False
        #check list 
        if person in one_gene:
            person_gene = 1
        elif person in two_genes: 
            person_gene = 2
       
            
     
        if person in have_trait:
            person_trait = True
        
        if people.get(person).get("mother") is None:
            work_dict_first.update({person : ( person_gene, person_trait, PROBS.get("gene")[person_gene] * PROBS.get("trait")[person_gene][person_trait]  )})
        else:
            pont = PROBS.get("trait")[person_gene][person_trait]
            work_dict_second.update({person : ( person_gene, person_trait,  pont )})


        
    for person in work_dict_second:
        get_gene = 0  #get gene probability
        if (people.get(person).get("mother") in work_dict_first.keys()) and (people.get(person).get("father") in work_dict_first.keys()):
            mother_genes = work_dict_first.get(people.get(person).get("mother"))[0]
            father_genes = work_dict_first.get(people.get(person).get("father"))[0]
            child_genes = work_dict_second.get(person)[0]
            if child_genes == 0:
                get_gene = transfer_table[mother_genes][0] * transfer_table[father_genes][0] 
            elif child_genes == 1:
                get_gene = transfer_table[mother_genes][1] * transfer_table[father_genes][0] + transfer_table[mother_genes][0] * transfer_table[father_genes][1]
            else: 
                get_gene = transfer_table[mother_genes][1] * transfer_table[father_genes][1] 
            
            work_dict_first.update({person : (child_genes, work_dict_second.get(person)[1] , work_dict_second.get(person)[2] * get_gene)})
            


    joint_prob = 1
    for person in work_dict_first:
        joint_prob = joint_prob * work_dict_first.get(person)[2]

    
       
    return joint_prob


    raise NotImplementedError


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """

    for person in probabilities:
        person_genes = 0
        person_trait = False
        
        if person in one_gene:
            person_genes = 1
        elif person in two_genes:
            person_genes = 2

        if person in have_trait:
            person_trait = True

          
        probabilities[person]["gene"][person_genes] = p + probabilities[person]["gene"][person_genes]
         
        probabilities[person]["trait"][person_trait] = p + probabilities[person]["trait"][person_trait]
       

        

    return
    raise NotImplementedError


def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities:
        #normalize gene distribution:
        val0 = probabilities[person]["gene"][0]
        val1 = probabilities[person]["gene"][1]
        val2 = probabilities[person]["gene"][2]

        probabilities[person]["gene"][0] = round(val0 / (val0 + val1 + val2),4)
        probabilities[person]["gene"][1] = round(val1 / (val0 + val1 + val2),4)
        probabilities[person]["gene"][2] = round(val2 / (val0 + val1 + val2),4)

        #normalize trait distribution:
        valTrue = probabilities[person]["trait"][True]
        valFalse = probabilities[person]["trait"][False]
        probabilities[person]["trait"][True] = round(valTrue / (valFalse + valTrue),4)
        probabilities[person]["trait"][False] = round(valFalse / (valFalse + valTrue),4)

    
    return
    raise NotImplementedError


if __name__ == "__main__":
    main()
