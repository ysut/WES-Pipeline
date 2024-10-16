/*
Copyright © 2024 NAME HERE <EMAIL ADDRESS>

*/
package cmd

import (
	"fmt"
	"github.com/c-bata/go-prompt"
	"github.com/spf13/cobra"
)

var date, name string

// probandCmd represents the proband command
var probandCmd = &cobra.Command{
	Use:   "proband",
	Short: "For proband-only pedigree format file.",
	Long: `For proband-only pedigree format file.
and usage of using your command. For example:

pedigo proband --output ./output.ped
`,

	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println("proband called")
	},
}

func init() {
	rootCmd.AddCommand(probandCmd)

	addCmd.Flags().StringVarP(&date, "date", "D", "", "Anniversary date in the format 'YYYYMMDD'")
	addCmd.Flags().StringVar(&name, "name", "", "Anniversary date name")


	// Here you will define your flags and configuration settings.

	// Cobra supports Persistent Flags which will work for this command
	// and all subcommands, e.g.:
	// probandCmd.PersistentFlags().String("foo", "", "A help for foo")

	// Cobra supports local flags which will only run when this command
	// is called directly, e.g.:
	// probandCmd.Flags().BoolP("toggle", "t", false, "Help message for toggle")
}
