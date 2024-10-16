package main

import (
	"fmt"
	"bufio"

	"github.com/c-bata/go-prompt"
)

func completer(d prompt.Document) []prompt.Suggest {
	s := []prompt.Suggest{
		{Text: "users", Description: "S"},
		{Text: "articles", Description: "Store the article text posted by user"},
		{Text: "comments", Description: "Store the text commented to articles"},
	}
	return prompt.FilterHasPrefix(s, d.GetWordBeforeCursor(), true)
}

func sexCompleter(d prompt.Document) []prompt.Suggest {
	s := []prompt.Suggest{
		{Text: "2", Description: "Female"},
		{Text: "1", Description: "Male"},
		{Text: "0", Description: "NA for sex or gender information."},
	}
	return prompt.FilterHasPrefix(s, d.GetWordBeforeCursor(), true)
}

func main() {
	fmt.Println("Please enter a Family ID.")
	familyID := prompt.Input("> ", completer)
	fmt.Println("Faimly ID: " + familyID)

	fmt.Println("Please enter a father ID.")
	fmt.Println("Please enter a mother ID.")
	fmt.Println("Please select sex.")
	form := prompt.Input("> ", sexCompleter)
	fmt.Println("Sex: " + form) // <- codeを性別にマッピングしたやつからとる

	
	// fmt.Println("Please select phenotype.")

}

// go build -o ./main cmd/main.go && ./main


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