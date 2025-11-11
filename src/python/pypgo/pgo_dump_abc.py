import pypgo
import sys

if __name__ == "__main__":
    print(pypgo.convert_animation_to_abc(sys.argv[1], sys.argv[2]))