import argparse
import csv
import glob
import os.path


METHODS = ("FACS", "droplet")
INPUT_ANNOTATIONS_DIR = "../03_tissue_annotation_csv"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", type=str, choices=METHODS, required=True)
    return parser.parse_args()


def get_header(csv_filename):
    with open(csv_filename) as csv_file:
        header = next(csv.reader(csv_file))
    return header


def main():
    args = parse_args()

    all_annotations_csv_name = f"annotations_{args.method}.csv"
    if os.path.isfile(all_annotations_csv_name):
        raise RuntimeError(f"{all_annotations_csv_name} already exists. Please remove before running.")
    input_annotations = glob.glob(f"{INPUT_ANNOTATIONS_DIR}/*_{args.method.lower()}_annotation.csv")

    # Get the header from one of the annotation files so we don't have to
    # worry about if we have processed the header or not in the loop.
    # This will also allow us to verify all of the headers are the same.
    first_header_seen = get_header(input_annotations[0])

    with open(all_annotations_csv_name, 'w') as all_annotations_csv_file:
        all_annotations_csv = csv.writer(all_annotations_csv_file, delimiter=',', quoting=csv.QUOTE_ALL, quotechar='"')
        all_annotations_csv.writerow(first_header_seen)

        for annotation_file in input_annotations:
            with open(annotation_file) as csv_file:
                annotation_data = csv.reader(csv_file, delimiter=',')
                annotation_header = next(annotation_data)
                if annotation_header != first_header_seen:
                    os.remove(all_annotations_csv_name)
                    raise RuntimeError(
                        f"""The annotation file: {annotation_file} does not match the first header seen: {first_header_seen}.
                           Failed to create {all_annotations_csv_name}.
                        """
                    )
                for row in annotation_data:
                    all_annotations_csv.writerow(row)


if __name__ == "__main__":
    main()
