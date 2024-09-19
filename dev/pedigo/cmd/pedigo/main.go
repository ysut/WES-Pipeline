package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"

	"github.com/olekukonko/tablewriter"
)

type Individual struct {
	FamilyID     string
	IndividualID string
	FatherID     string
	MotherID     string
	Sex          string
	Phenotype    string
}

func exists(name string) bool {
	_, err := os.Stat(name)
	// If the file does not exist, os.IsNotExist(err) returns true
	return !os.IsNotExist(err)
}

func main() {
	// Flags
	var outputFlag = flag.String("output", "./output.ped", "output file path")
	var yesFlag = flag.Bool("yes", false, "Non-confirmation mode")
	flag.Parse()

	scanner := bufio.NewScanner(os.Stdin)
	// Check if the output file exists
	if exists(*outputFlag) {
		fmt.Println("The output file already exists. Do you want to overwrite it? [Y/n]")
		var inOverwrite string
		for {
			scanner.Scan()
			inOverwrite = scanner.Text()
			if inOverwrite == "Y" || inOverwrite == "n" {
				break
			}
			fmt.Println("Please select one of [Y/n]")
		}
		if inOverwrite == "n" {
			log.Fatal("Please provide a different output file path using --output flag")
		} else {
			// Remove the existing file
			err := os.Remove(*outputFlag)
			if err != nil {
				log.Fatal(err)
			}
		}
	}

	// Create a new file
	file, err := os.Create(*outputFlag)
	if err != nil {
		log.Println(err)
	}
	defer file.Close()

	n := checkTotalFamilies(scanner)

	// Iterate through each project and add individuals to the array
	var individuals []Individual
	for i := 0; i < n; i++ {
		fmt.Printf("\n== Please provide the information for family #%d (total %d) ==", i+1, n)
		// STEP 1 - Get the family ID
		familyID := enterFamilyID(scanner)

		// STEP2 - Get the number of individuals in the family
		mode, numIndivisuals := checkNumFamilyMembers(scanner)
		var isConfirmed bool = false
		switch mode {
		case "p":
			for !isConfirmed {
				fmt.Println("\nProvide the proband information")
				_, _ = enterSingleIndividualInfo(scanner, familyID, &individuals, mode)
				isConfirmed = confirmCurrentInputs(scanner, &individuals, familyID, yesFlag)
			}
			continue
		case "t":
			for !isConfirmed {
				fmt.Println("\nProvide the proband information")
				fatherID, motherID := enterSingleIndividualInfo(scanner, familyID, &individuals, mode)
				enterParentsInfo(scanner, familyID, &individuals, fatherID, motherID)
				isConfirmed = confirmCurrentInputs(scanner, &individuals, familyID, yesFlag)
			}
			continue
		case "q":
			for !isConfirmed {
				fmt.Println("\nProvide the proband information")
				fatherID, motherID := enterSingleIndividualInfo(scanner, familyID, &individuals, mode)
				enterParentsInfo(scanner, familyID, &individuals, fatherID, motherID)
				enterSiblingInfo(scanner, familyID, &individuals, fatherID, motherID)
				isConfirmed = confirmCurrentInputs(scanner, &individuals, familyID, yesFlag)
			}
			continue
		case "n":
			for !isConfirmed {
				for j := 0; j < numIndivisuals; j++ {
					fmt.Printf("\nProvide the information for individual #%d (total %d)", j+1, numIndivisuals)
					_, _ = enterSingleIndividualInfo(scanner, familyID, &individuals, mode)
				}
				isConfirmed = confirmCurrentInputs(scanner, &individuals, familyID, yesFlag)
			}
			continue
		}
	}

	// Write to the output file

	displayAllIndividuals(&individuals)

	log.Println("Saving a pedigree foramt file...")
	for _, individual := range individuals {
		line := fmt.Sprintf("%s\t%s\t%s\t%s\t%s\t%s\n",
			individual.FamilyID,
			individual.IndividualID,
			individual.FatherID,
			individual.MotherID,
			individual.Sex,
			individual.Phenotype)
		file.WriteString(line)
	}
	log.Printf("The file has been saved successfully\n")
}

//-------------------//
//     Functions     //
//-------------------//

func checkTotalFamilies(scanner *bufio.Scanner) int {
	// var inTotalFamilies string
	var n int
	for {
		fmt.Println("Enter a total number of families: ")
		scanner.Scan()
		inTotalFamilies := scanner.Text()
		n, _ = strconv.Atoi(inTotalFamilies)
		if n > 0 {
			break
		}
		fmt.Println("Please enter a natural number")
	}
	return n
}

func enterFamilyID(scanner *bufio.Scanner) string {
	var familyID string
	for {
		fmt.Println("\nEnter the family ID: ")
		scanner.Scan()
		familyID = scanner.Text()
		if len(familyID) > 0 {
			break
		}
	}
	return familyID
}

func checkNumFamilyMembers(scanner *bufio.Scanner) (string, int) {
	var mode string
	var numIndivisuals int

	for {
		fmt.Println("\nHow many individuals are in this family?: ")
		fmt.Println("(Shortcut) p for proband, t for trio, and q for quad")

		scanner.Scan()
		inStringIndividuals := scanner.Text()

		// When the input is a shortcut
		switch inStringIndividuals {
		case "p":
			mode = inStringIndividuals
			numIndivisuals = 1
			goto end
		case "t":
			mode = inStringIndividuals
			numIndivisuals = 3
			goto end
		case "q":
			mode = inStringIndividuals
			numIndivisuals = 4
			goto end
		default:
			goto L
		}
	L: // When the input is a natural number
		numIndivisuals, _ = strconv.Atoi(inStringIndividuals)
		if numIndivisuals > 0 {
			mode = "n"
			goto end
		}
		fmt.Println("Please select one of [a natural number, p , t, q]")
	}
end:
	return mode, numIndivisuals
}

func enterIndividualID(scanner *bufio.Scanner) string {
	scanner.Scan()
	individualID := scanner.Text()
	return individualID
}

func enterSexInfo(scanner *bufio.Scanner) string {
	var inSexInfo string
	var n int
	for {
		scanner.Scan()
		inSexInfo = scanner.Text()
		n, _ = strconv.Atoi(inSexInfo)
		if n < 3 {
			break
		}
		fmt.Println("Please select one of [1=male; 2=female; other=unknown]")
	}
	return inSexInfo
}

func enterPhenotypeInfo(scanner *bufio.Scanner) string {
	var inPhenotypeInfo string
	var n int
	for {
		scanner.Scan()
		inPhenotypeInfo = scanner.Text()
		n, _ = strconv.Atoi(inPhenotypeInfo)
		if n == -9 || n == 0 || n == 1 || n == 2 {
			break
		}
		fmt.Println("Please select one of [1=male; 2=female; other=unknown]")
		fmt.Println("(Optional) You can provide not listed above but HPO IDs")
	}
	return inPhenotypeInfo
}

func enterSingleIndividualInfo(
	scanner *bufio.Scanner, familyID string, individuals *[]Individual, mode string) (string, string) {

	fmt.Println("Enter the individual ID: ")
	probandID := enterIndividualID(scanner)

	var fatherID string
	var motherID string
	if mode == "p" {
		fatherID = "0"
		motherID = "0"
	} else {
		fmt.Println("\nEnter the Father ID: ")
		fatherID = enterIndividualID(scanner)
		fmt.Println("\nEnter the Mother ID: ")
		motherID = enterIndividualID(scanner)
	}
	fmt.Println("\nSelect one of [1=male; 2=female; other=unknown]")
	form := enterSexInfo(scanner)
	fmt.Println("\nEnter one of [-9=missing; 0=missing; 1=unaffected; 2=affected]")
	fmt.Println("(Optional) You can also provide HPO IDs not listed above")
	phenotype := enterPhenotypeInfo(scanner)

	individual := Individual{
		FamilyID:     familyID,
		IndividualID: probandID,
		FatherID:     fatherID,
		MotherID:     motherID,
		Sex:          form,
		Phenotype:    phenotype,
	}
	*individuals = append(*individuals, individual)

	return fatherID, motherID
}

func enterParentsInfo(
	scanner *bufio.Scanner, familyID string, individuals *[]Individual, fatherID string, motherID string) {

	fmt.Println("\nEnter the father's phenotype")
	fmt.Println("Enter one of [-9=missing; 0=missing; 1=unaffected; 2=affected]")
	fmt.Println("(Optional) You can also provide HPO IDs not listed above")
	dadPhenotype := enterPhenotypeInfo(scanner)

	father := Individual{
		FamilyID:     familyID,
		IndividualID: fatherID,
		FatherID:     "0",
		MotherID:     "0",
		Sex:          "1",
		Phenotype:    dadPhenotype,
	}
	*individuals = append(*individuals, father)

	fmt.Println("\nEnter the mother's phenotype")
	fmt.Println("Enter one of [-9=missing; 0=missing; 1=unaffected; 2=affected]")
	fmt.Println("(Optional) You can also provide HPO IDs not listed above")
	momPhenotype := enterPhenotypeInfo(scanner)

	mother := Individual{
		FamilyID:     familyID,
		IndividualID: motherID,
		FatherID:     "0",
		MotherID:     "0",
		Sex:          "2",
		Phenotype:    momPhenotype,
	}
	*individuals = append(*individuals, mother)
}

func enterSiblingInfo(
	scanner *bufio.Scanner, familyID string, individuals *[]Individual, fatherID string, motherID string) {

	fmt.Println("\nEnter the sibling information (ID, Sex, and Phenotype)")
	fmt.Println("Enter the Sibling ID: ")
	probandID := enterIndividualID(scanner)

	fmt.Println("\nSelect one of [1=male; 2=female; other=unknown]")
	form := enterSexInfo(scanner)

	fmt.Println("\nEnter one of [-9=missing; 0=missing; 1=unaffected; 2=affected]")
	fmt.Println("(Optional) You can also provide HPO IDs not listed above")
	sibPhenotype := enterPhenotypeInfo(scanner)

	sibling := Individual{
		FamilyID:     familyID,
		IndividualID: probandID,
		FatherID:     fatherID,
		MotherID:     motherID,
		Sex:          form,
		Phenotype:    sibPhenotype,
	}
	*individuals = append(*individuals, sibling)
}

func displayCurrentIndividuals(individuals *[]Individual, familyID string) [][]string {
	var displayData [][]string
	for _, individual := range *individuals {
		if individual.FamilyID == familyID {
			displayData = append(displayData, []string{individual.FamilyID, individual.IndividualID, individual.FatherID, individual.MotherID, individual.Sex, individual.Phenotype})
		}
	}
	table := tablewriter.NewWriter(os.Stdout)
	table.SetHeader([]string{"Family ID", "Individual ID", "Father ID", "Mother ID", "Sex", "Phenotype"})
	for _, v := range displayData {
		table.Append(v)
	}
	table.Render()

	return displayData
}

func displayAllIndividuals(individuals *[]Individual) {
	fmt.Println("\n>>> All family data are shown as below")
	var displayData [][]string
	for _, individual := range *individuals {
		displayData = append(displayData, []string{individual.FamilyID, individual.IndividualID, individual.FatherID, individual.MotherID, individual.Sex, individual.Phenotype})
	}
	table := tablewriter.NewWriter(os.Stdout)
	table.SetHeader([]string{"Family ID", "Individual ID", "Father ID", "Mother ID", "Sex", "Phenotype"})
	table.SetCaption(true, "\n")
	for _, v := range displayData {
		table.Append(v)
	}
	table.Render()
}

func confirmCurrentInputs(scanner *bufio.Scanner, individuals *[]Individual, familyID string, yesFlag *bool) bool {
	// Confirm the current inputs and decide whether to write to the output file
	if *yesFlag {
		return true
	} else {
		displayData := displayCurrentIndividuals(individuals, familyID)
		fmt.Println("Do you want to write above information to the output file? [Y/n]")
		var inYn string
		for {
			scanner.Scan()
			inYn = scanner.Text()
			if inYn == "Y" || inYn == "y" || inYn == "N" || inYn == "n" {
				break
			}
		}
		if inYn == "n" || inYn == "N" {
			// Remove the last family's individuals and return false
			tmpArray := []Individual{}
			*individuals = append(tmpArray, (*individuals)[:len(*individuals)-len(displayData)]...)
			return false
		} else {
			return true
		}
	}
}
